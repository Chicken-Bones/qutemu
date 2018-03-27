#ifndef CHASTE_SLABSIMULATION_H
#define CHASTE_SLABSIMULATION_H

#include <MonodomainProblem.hpp>
#include <QuadraticMesh.hpp>
#include <SimpleStimulus.hpp>
//#include <ORdCvode.hpp>
#include <ORd2011endo_fkatpCvode.hpp>
//#include <BeelerReuter1977Cvode.hpp>

class CornerStimulus2dCellFactory : public AbstractCardiacCellFactory<2>
{
private:
    double corner_size;
    boost::shared_ptr<SimpleStimulus> mpStimulus;

public:
    CornerStimulus2dCellFactory(double corner_size, double amp, double dur)
            : AbstractCardiacCellFactory<2>(),
              corner_size(corner_size),
              mpStimulus(new SimpleStimulus(amp, dur))
    {
    }

    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<2>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        boost::shared_ptr<AbstractStimulusFunction> stim =
                x<corner_size+1e-6 && y<corner_size+1e-6 ? //assumes x and y > 0
                (boost::shared_ptr<AbstractStimulusFunction>) mpStimulus : mpZeroStimulus;

        //AbstractCvodeCell* p_cell = new CellBeelerReuter1977FromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver>(), stim);
        //AbstractCvodeCell* p_cell = new CellORdFromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver>(), stim);
        AbstractCvodeCell* p_cell = new CellORd2011endo_fkatpFromCellMLCvode(boost::shared_ptr<AbstractIvpOdeSolver>(), stim);
        p_cell->SetTolerances(1e-5, 1e-7);
        return p_cell;
    }
};

class SlabSimulation
{
public:
    static void BasicSquare(int n = 400, double size = 10, double duration_ms = 500, int timestep_us = 100) {
        HeartConfig::Instance()->SetSimulationDuration(duration_ms);
        std::stringstream ss;
        ss << "qutemu/BasicMonodomainMesh/" << n << "n" << size << "mm@" << timestep_us << "us";
        HeartConfig::Instance()->SetOutputDirectory(ss.str());
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        CornerStimulus2dCellFactory cell_factory(0.5, -60.0 * 2000.0, 5.0);
        MonodomainProblem<2> problem( &cell_factory );

        AbstractTetrahedralMesh<2, 2> *mesh;
        if (PetscTools::IsParallel())
            mesh = new DistributedTetrahedralMesh<2, 2>();
        else
            mesh = new TetrahedralMesh<2, 2>();

        mesh->ConstructRegularSlabMesh(size/n, size, size);
        problem.SetMesh(mesh);

        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000.0);
        HeartConfig::Instance()->SetCapacitance(1.0);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.5, 1.5));
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(timestep_us/1000.0, timestep_us/1000.0, 1.0);

        problem.Initialise();
        problem.Solve();

        delete mesh;
    }
};

#endif //CHASTE_SLABSIMULATION_HPP
