#include "ExecutableSupport.hpp"
#include "CommandLineArguments.hpp"
#include "MonodomainProblem.hpp"
#include "SimpleStimulus.hpp"
#include "SteadyStateRunner.hpp"
#include "CardiacSimulationArchiver.hpp"

#include "Maleckar2008_LA_1h2HzCvodeOpt.hpp"
#include "Maleckar2008_RA_1h2HzCvodeOpt.hpp"
#include "courtemanche_ramirez_nattel_1998_SRCvodeOpt.hpp"
#include "courtemanche_ramirez_nattel_1998_cAFCvodeOpt.hpp"

#include "QutemuLog.hpp"
#include "QutemuVersion.hpp"
#include "ConductivityReader.hpp"
#include "ActivationMapOutputModifier.hpp"
#include "TimedStimulus.hpp"

#include <sys/resource.h>
#include <Version.hpp>
#include <boost/lexical_cast.hpp>

class AtrialCellFactory : public AbstractCardiacCellFactory<3> // <3> here
{
public:
    static const int MALECKAR = 0;
    static const int COURTEMANCHE_SR = 1;
    static const int COURTEMANCHE_CAF = 2;

private:

    boost::shared_ptr<AbstractStimulusFunction> p_stim_sinus;
    boost::shared_ptr<AbstractStimulusFunction> p_stim_extra;
    int p_cell_model;

public:
    AtrialCellFactory() {}

    AtrialCellFactory(boost::shared_ptr<AbstractStimulusFunction> p_stim_sinus, boost::shared_ptr<AbstractStimulusFunction> p_stim_extra, int p_cell_model) :
            AbstractCardiacCellFactory<3>(),
	        p_stim_sinus(p_stim_sinus),
	        p_stim_extra(p_stim_extra),
            p_cell_model(p_cell_model)
    {
        if (p_cell_model < MALECKAR || p_cell_model > COURTEMANCHE_CAF)
            EXCEPTION("Unknown Cell Model " << p_cell_model);
    }
    
    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        unsigned lvrv = 1;
        unsigned pacing_site = 0;
        if (pNode->HasNodeAttributes()) {//halo nodes have no attributes
            std::vector<double>& attributes = pNode->rGetNodeAttributes();
            lvrv = (unsigned)attributes[0];
            pacing_site = (unsigned)attributes[1];
        }

        if (lvrv < 1 || lvrv > 2)
            EXCEPTION("invalid lvrv " << lvrv << " at node " << pNode->GetIndex());

        boost::shared_ptr<AbstractStimulusFunction> stimulus;
        switch (pacing_site) {
            case 1:
                stimulus = p_stim_sinus;
                break;
            case 2:
                stimulus = p_stim_extra;
                break;
            case 0:
                stimulus = mpZeroStimulus;
                break;
            default:
                EXCEPTION("Unknown Pacing Site " << pacing_site << " at node " << pNode->GetIndex());
        }

        switch (p_cell_model) {
            case MALECKAR:
                if (lvrv == 1)
                    return new CellMaleckar2008_LA_1h2HzFromCellMLCvodeOpt(mpSolver, stimulus);
                else
                    return new CellMaleckar2008_RA_1h2HzFromCellMLCvodeOpt(mpSolver, stimulus);
            case COURTEMANCHE_SR:
                return new Cellcourtemanche_ramirez_nattel_1998_SRFromCellMLCvodeOpt(mpSolver, stimulus);
            case COURTEMANCHE_CAF:
                return new Cellcourtemanche_ramirez_nattel_1998_cAFFromCellMLCvodeOpt(mpSolver, stimulus);
            default:
                EXCEPTION("Um");

        }
    }
};


class AtrialConductivityModifier : public AbstractConductivityModifier<3,3>
{
private:
	
    c_matrix<double,3,3> mTensor;
    AbstractTetrahedralMesh<3,3>* pMesh;
    std::vector<float> conductivities;
    
public:
    AtrialConductivityModifier() {}

    AtrialConductivityModifier(AbstractTetrahedralMesh<3,3>* mesh, const std::vector<float> &conductivities) :
            AbstractConductivityModifier<3,3>(),
            mTensor(zero_matrix<double>(3,3)),
            pMesh(mesh),
            conductivities(conductivities)
    {
        assert(conductivities.size() == 0 || conductivities.size() == mesh->GetNumElements());
    }
    
    c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,3,3>& rOriginalConductivity, unsigned domainIndex)
    {
        double gll;
        double gtt;

    	// read element attribute and tune conductivity
    	int a=pMesh->GetElement(elementIndex)->rGetElementAttributes()[0];


    	//*************** TUNE CONDUCTIVITY PER TISSUE CLASS !! ******************/
    	if ((a==32) || (a==33) || (a==76) || (a==111) || (a==192) || (a==193) || (a==194) || (a==195) || (a==196) || (a==197) || (a==198) || (a==199) || (a==112) || (a==161) || (a==162))       
    	{
    		gll = 2.0;
    		gtt = 1.0;
    	}
    	// sinus node and surroundings
    	else if ((a==80) || (a==81) || (a==82) || (a==83) || (a==84) || (a==85) || (a==86))
    	{
    		gll = 0.5;
    		gtt = 0.5;
    	}
    	
    	else if (a==72) // crista terminalis
    	{
    		gll = 7.7;
    		gtt = 0.7;
    	}
    	
    	else if ((a==74) || (a==98) || (a==102) || (a==103)) // pectinatae muscles, Bachman Bundle
    	{
    		gll = 5.5;
    		gtt = 2.75;
    	}
    	else if (a==104) // right atrial ithmus
    	{
    		gll = 1.1;
    		gtt = 1.1;
    	}
        else
        {
            EXCEPTION("Unknown cell class " << a << " at " << elementIndex);
        }

        if (conductivities.size() > 0)
        {
            double f = conductivities[elementIndex];
            gll *= f;
            gtt *= f;
        }

    	mTensor(0,0) = gll;
        mTensor(1,1) = gtt;
        mTensor(2,2) = gtt;
        return mTensor;
    }
};

class AtrialFibrosis
{
private:
    double GetMemoryUsage()
    {
    	struct rusage rusage;
    	getrusage( RUSAGE_SELF, &rusage );
    	return (double)(rusage.ru_maxrss)/(1024);// Convert KB to MB
    }

    void SetSchemaLocations()
    {
        std::string root_dir = boost::filesystem::current_path().append("xsd/", boost::filesystem::path::codecvt()).string();
        std::map<std::string, std::string> schemaLocations;
        schemaLocations[""] = root_dir + "ChasteParameters_1_1.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_0"] = root_dir + "ChasteParameters_2_0.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_1"] = root_dir + "ChasteParameters_2_1.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_2"] = root_dir + "ChasteParameters_2_2.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/2_3"] = root_dir + "ChasteParameters_2_3.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_0"] = root_dir + "ChasteParameters_3_0.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_1"] = root_dir + "ChasteParameters_3_1.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_3"] = root_dir + "ChasteParameters_3_3.xsd";
        schemaLocations["https://chaste.comlab.ox.ac.uk/nss/parameters/3_4"] = root_dir + "ChasteParameters_3_4.xsd";
        HeartConfig::Instance()->SetFixedSchemaLocations(schemaLocations);
    }

    void OverrideVoltageLookupRange() {
        boost::shared_ptr<AbstractIvpOdeSolver> noSolver;
        boost::shared_ptr<AbstractStimulusFunction> noStim;
        AbstractLookupTableCollection *tables[] = {
                CellMaleckar2008_LA_1h2HzFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                CellMaleckar2008_RA_1h2HzFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                Cellcourtemanche_ramirez_nattel_1998_SRFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                Cellcourtemanche_ramirez_nattel_1998_cAFFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
        };

        for (int i = 0; i < 4; i++) {
            AbstractLookupTableCollection *p_tables = tables[i];
            double min, step, max;
            p_tables->GetTableProperties("membrane__V", min, step, max);
            min = -150.0001;
            p_tables->SetTableProperties("membrane__V", min, step, max);
            p_tables->RegenerateTables();
        }
    }

    template<class T>
    std::vector<T> ParseMultiValueOption(std::string optname)
    {
        std::string option = CommandLineArguments::Instance()->GetValueCorrespondingToOption(optname);
        std::vector<T> nodes;
        if (option.find(',') != std::string::npos) {
            std::stringstream ss(option);
            std::string item;
            try {
                while (std::getline(ss, item, ','))
                    nodes.push_back(boost::lexical_cast<T>(item));
            } catch (boost::bad_lexical_cast e) {
                EXCEPTION("Invalid value: '" << item << "' while parsing " << optname);
            }

            return nodes;
        }

        try
        {
            nodes.push_back(boost::lexical_cast<T>(option));
            return nodes;
        }
        catch(const boost::bad_lexical_cast &e)
        {}

        std::ifstream file(option.c_str(), std::ios::in);
        if (!file.is_open())
            EXCEPTION("Couldn't open file: " + option);

        std::string line;
        try {
            while (std::getline(file, line))
                nodes.push_back(boost::lexical_cast<T>(line));
        } catch (boost::bad_lexical_cast e) {
            EXCEPTION("Invalid value: '" << line << "' in file '" << option << "' while parsing " << optname);
        }

        return nodes;
    }

    double GetDoubleOption(std::string pname, double defaultValue) {
        return CommandLineArguments::Instance()->OptionExists(pname) &&
               CommandLineArguments::Instance()->GetNumberOfArgumentsForOption(pname) > 0 ?
               CommandLineArguments::Instance()->GetDoubleCorrespondingToOption(pname) :
               defaultValue;
    }

    int GetIntOption(std::string pname, int defaultValue) {
        return CommandLineArguments::Instance()->OptionExists(pname) &&
               CommandLineArguments::Instance()->GetNumberOfArgumentsForOption(pname) > 0 ?
               CommandLineArguments::Instance()->GetIntCorrespondingToOption(pname) :
               defaultValue;
    }

    OutputFileHandler InitOutput() {
        if (CommandLineArguments::Instance()->OptionExists("-outdir"))
            HeartConfig::Instance()->SetOutputDirectory(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-outdir"));

        LOG("outdir: " << HeartConfig::Instance()->GetOutputDirectory());
        return OutputFileHandler(HeartConfig::Instance()->GetOutputDirectory(), false);
    }

    void InitTimesteps() {
        double duration = GetDoubleOption("-duration", 5.0);
        HeartConfig::Instance()->SetSimulationDuration(duration); //ms

        double odet = GetDoubleOption("-odet", 0.02);
        double pdet = GetDoubleOption("-pdet", odet);
        double interval = GetDoubleOption("-interval", 5.0);
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(odet, pdet, interval);

        LOG("** TIMESTEPS **");
        LOG("duration: " << duration << "ms");
        LOG("interval: " << interval << "ms");
        LOG("odet    : " << odet << "ms");
        LOG("pdet    : " << pdet << "ms");
    }

    boost::shared_ptr<AbstractStimulusFunction> InitStimulus(std::string name, std::vector<double> &rStimTimes,
            unsigned defaultN, double defaultDelay, double defaultPeriod) {

        bool user_defined = CommandLineArguments::Instance()->OptionExists("-"+name);

        double offset = rStimTimes.empty() ? 0 : rStimTimes.back();
        offset += GetDoubleOption("-d"+name, user_defined ? 0 : defaultDelay);

        double period;
        std::vector<double> times;
        if (user_defined) {
            times = ParseMultiValueOption<double>("-"+name);
            for (double &t : times)
                t += offset;
        } else {
            period = GetDoubleOption("-p"+name, defaultPeriod);
            int nstimuli = (unsigned)GetIntOption("-n"+name, defaultN);
            for (int i = 0; i < nstimuli; i++)
                times.push_back(offset + i*period);
        }

        LOG(name << ":");
        LOG("\tcycles: " << times.size());
        if (!times.empty()) {
            LOG("\tfirst : " << times[0] << "ms");
            LOG("\tlast  : " << times.back() << "ms");
            if (user_defined) {
                std::stringstream ss;
                ss << times[0];
                for (unsigned i = 1; i < times.size(); i++)
                    ss << ", " << times[i];

                LOG("\ttimes : " << ss.str());
            }
            else
                LOG("\tperiod: " << period << "ms");
        }

        rStimTimes.insert( rStimTimes.end(), times.begin(), times.end() );
        return boost::shared_ptr<AbstractStimulusFunction>(new TimedStimulus(-80000.0, 1.0, times));
    }

    AtrialCellFactory InitCellFactory(std::vector<double> &rStimTimes) {
        LOG("** CELLS **")
        std::string cellopt = "maleckar";
        if (CommandLineArguments::Instance()->OptionExists("-cell"))
            cellopt = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-cell");

        int cell_model;
        if (cellopt == "maleckar")
            cell_model = AtrialCellFactory::MALECKAR;
        else if (cellopt == "courtemanche_sr")
            cell_model = AtrialCellFactory::COURTEMANCHE_SR;
        else if (cellopt == "courtemanche_caf")
            cell_model = AtrialCellFactory::COURTEMANCHE_CAF;
        else
            EXCEPTION("Unknown Cell Model: " << cellopt);
        LOG("cell: " << cellopt);

        auto p_stim_sinus = InitStimulus("sinus", rStimTimes, 4, 0, 500);
        auto p_stim_extra = InitStimulus("extra", rStimTimes, 6, 400, 300);

        OverrideVoltageLookupRange();
        return AtrialCellFactory(p_stim_sinus, p_stim_extra, cell_model);
    }

    void ApplyPerm(std::vector<unsigned int> &nodes, const std::vector<unsigned int> &permutation) {
        if (permutation.size() > 0) {
            for (unsigned i = 0; i < nodes.size(); i++)
                nodes[i] = permutation[nodes[i]];
        }

        std::sort(nodes.begin(), nodes.end());
    }

    void AddActivationMap(MonodomainProblem<3> *problem, const std::vector<double> &rStimTimes) {
        if (CommandLineArguments::Instance()->OptionExists("-nosnapshots"))
            return;

        //get initial value from cell system
        unsigned local_node0 = problem->rGetMesh().GetDistributedVectorFactory()->GetLow();
        AbstractCardiacCellInterface* cell = problem->GetMonodomainTissue()->GetCardiacCell(local_node0);
        auto* system = dynamic_cast<AbstractUntemplatedParameterisedSystem *>(cell);
        double resting = system->GetSystemInformation()->GetInitialConditions()[cell->GetVoltageIndex()];
        double threshold = GetDoubleOption("-activation", -40);

        problem->AddOutputModifier(boost::shared_ptr<AbstractOutputModifier>(new ActivationMapOutputModifier("snapshots.h5", threshold, resting, rStimTimes)));
        problem->AddOutputModifier(boost::shared_ptr<AbstractOutputModifier>(new ActivationMapOutputModifier("snapshots_dyn.h5", threshold, resting)));

        LOG("activationmap:")
        LOG("\tthreshold: " << threshold << "mV")
        LOG("\tresting  : " << resting << "mV")
    }

    MonodomainProblem<3> *InitProblem(AtrialCellFactory *cell_factory)
    {
        CommandLineArguments* args = CommandLineArguments::Instance();
        HeartConfig* heartConfig = HeartConfig::Instance();
        MonodomainProblem<3>* problem;

        if (args->OptionExists("-svi"))
            heartConfig->SetUseStateVariableInterpolation(true);

        LOG("** PROBLEM **")
        if (args->OptionExists("-loaddir")) {
            std::string loaddir = args->GetStringCorrespondingToOption("-loaddir");
            LOG("loaddir: " << loaddir);
            problem = CardiacSimulationArchiver<MonodomainProblem<3> >::Load(loaddir);
        }
        else {
            std::string meshfile = args->GetStringCorrespondingToOption("-meshfile");
            LOG("meshfile: " << meshfile);
            heartConfig->SetMeshFileName(meshfile, cp::media_type::Axisymmetric);

            problem = new MonodomainProblem<3>(cell_factory);
            problem->SetWriteInfo();
            problem->Initialise();
        }

        LOG("svi: " << (heartConfig->GetUseStateVariableInterpolation() ? "true" : "false"))

        if (args->OptionExists("-nodes")) {
            std::vector<unsigned> nodes = ParseMultiValueOption<unsigned>("-nodes");
            ApplyPerm(nodes, problem->rGetMesh().rGetNodePermutation());
            problem->SetOutputNodes(nodes);
            LOG("nodes: " << args->GetStringCorrespondingToOption("-nodes"));
        }

        heartConfig->SetOutputFilenamePrefix("results");
        heartConfig->SetVisualizeWithMeshalyzer(false);
        heartConfig->SetVisualizeWithCmgui(false);
        heartConfig->SetVisualizeWithVtk(false);
        heartConfig->SetVisualizeWithParallelVtk(args->OptionExists("-vtk"));
        LOG("vtk: " << (heartConfig->GetVisualizeWithParallelVtk() ? "true" : "false"));

        return problem;
    }

    void InitConductivities(MonodomainProblem<3> *problem, AtrialConductivityModifier& modifier) {
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(1.75, 0.19, 0.19));

        std::vector<float> conductivities;
        if (CommandLineArguments::Instance()->OptionExists("-condmod")) {
            std::string path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-condmod");
            LOG("conductivities: " << path);
            conductivities = ConductivityReader::ReadConductivities(FileFinder(path, RelativeTo::AbsoluteOrCwd));
        }

        modifier = AtrialConductivityModifier(&problem->rGetMesh(), conductivities);
        problem->GetTissue()->SetConductivityModifier( &modifier );
    }

    void Save(MonodomainProblem<3> *problem) {
        if (CommandLineArguments::Instance()->OptionExists("-savedir"))
        {
            std::string savedir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-savedir");
            LOG("savedir: " << savedir);
            CardiacSimulationArchiver<MonodomainProblem<3> >::Save(*problem, savedir);
        }
    }

    void WriteLog(OutputFileHandler out_dir)
    {
        if (!PetscTools::AmMaster())
            return;

        LOG("finished: " << std::setprecision(3) << std::fixed << Timer::GetWallTime() << "s");
        out_stream os = out_dir.OpenOutputFile("log.txt");
        *os << QutemuLog::GetLog();
        os->close();
    }

    void WritePermutation(OutputFileHandler out_dir, MonodomainProblem<3> *problem)
    {
        if (!PetscTools::AmMaster())
            return;

        std::vector<unsigned> perm_vec = problem->rGetMesh().rGetNodePermutation();
        if (perm_vec.size() == 0)
            return;

        out_stream os = out_dir.OpenOutputFile("permutation.txt");
        (*os) << "Meshfile Chaste" << std::endl;
        for (unsigned i = 0; i < perm_vec.size(); ++i)
            (*os) << i << ' ' << perm_vec[i] << std::endl;
    }

public:
    void TestMonodomain3dAtria() throw(Exception)
    {
        SetSchemaLocations();

        COUT("Initializing");
        OutputFileHandler out_dir = InitOutput();
        InitTimesteps();
        std::vector<double> stim_times;
        AtrialCellFactory cell_factory = InitCellFactory(stim_times);
        MonodomainProblem<3>* problem = InitProblem(&cell_factory);
        AtrialConductivityModifier conductivity_modifier;
        InitConductivities(problem, conductivity_modifier);
        AddActivationMap(problem, stim_times);

        COUT("Solving");
        problem->Solve();
        Save(problem);

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        if (PetscTools::AmMaster())
            out_dir.FindFile("progress_status.txt").Remove();

        WritePermutation(out_dir, problem);
        WriteLog(out_dir);

        delete problem;
        COUT("Success");
    }
};

int main(int argc, char *argv[])
{
    ExecutableSupport::InitializePetsc(&argc, &argv);
    ExecutableSupport::ShowParallelLaunching();

    if (PetscTools::AmMaster()) {
        std::cout << "This version of Chaste was compiled on:\n" << ChasteBuildInfo::GetBuildTime()
                  << " with build type " << ChasteBuildInfo::GetBuildInformation() << ".\n\n";

        std::cout << "This version of qutemu was compiled on:\n" << QutemuVersion::GetBuildTime()
                  << " with build type " << ChasteBuildInfo::GetBuildInformation() << ".\n\n";
    }

    int exit_code = ExecutableSupport::EXIT_OK;

    try
    {
        AtrialFibrosis().TestMonodomain3dAtria();
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
