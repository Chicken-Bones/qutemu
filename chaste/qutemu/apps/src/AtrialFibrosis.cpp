#include "ExecutableSupport.hpp"
#include "CommandLineArguments.hpp"
#include "MonodomainProblem.hpp"
#include "SimpleStimulus.hpp"
#include "SteadyStateRunner.hpp"
#include "CardiacSimulationArchiver.hpp"

#include "Maleckar2008_baseCvodeOpt.hpp"
#include "Maleckar2008_cAFCvodeOpt.hpp"
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

enum CellModel
{
    MALECKAR,
    MALECKAR_CAF,
    MALECKAR_ANNA,
    COURTEMANCHE_SR,
    COURTEMANCHE_CAF
};

template<unsigned DIM>
class AtrialCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    boost::shared_ptr<AbstractStimulusFunction> p_stim_sinus;
    boost::shared_ptr<AbstractStimulusFunction> p_stim_extra;
    int p_cell_model;

    using AbstractCardiacCellFactory<DIM>::mpSolver;
    using AbstractCardiacCellFactory<DIM>::mpZeroStimulus;

public:
    AtrialCellFactory() {}

    AtrialCellFactory(boost::shared_ptr<AbstractStimulusFunction> p_stim_sinus, boost::shared_ptr<AbstractStimulusFunction> p_stim_extra, int p_cell_model) :
            AbstractCardiacCellFactory<DIM>(),
            p_stim_sinus(p_stim_sinus),
            p_stim_extra(p_stim_extra),
            p_cell_model(p_cell_model)
    {
        if (p_cell_model < MALECKAR || p_cell_model > COURTEMANCHE_CAF)
            EXCEPTION("Unknown Cell Model " << p_cell_model);
    }
    
    AbstractCvodeCell* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
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
                return new CellMaleckar2008_baseFromCellMLCvodeOpt(mpSolver, stimulus);
            case MALECKAR_CAF:
                return new CellMaleckar2008_cAFFromCellMLCvodeOpt(mpSolver, stimulus);
            case MALECKAR_ANNA:
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

template<unsigned DIM>
class AtrialConductivityModifier : public AbstractConductivityModifier<DIM,DIM>
{
private:
	
    c_matrix<double,DIM,DIM> mTensor;
    AbstractTetrahedralMesh<DIM,DIM>* pMesh;
    std::vector<float> conductivities;
    
public:
    AtrialConductivityModifier() {}

    AtrialConductivityModifier(const std::vector<float> &conductivities) :
            AbstractConductivityModifier<DIM,DIM>(),
            mTensor(zero_matrix<double>(DIM,DIM)),
            pMesh(NULL),
            conductivities(conductivities)
    {
    }

    void SetMesh(AbstractTetrahedralMesh<DIM,DIM>* mesh) {
        pMesh = mesh;
        assert(conductivities.size() == 0 || conductivities.size() == mesh->GetNumElements());
    }

    void ApplyTissueConductivity(unsigned elementIndex, unsigned tissue_class) {
        double gll;
        double gtt;
        switch (tissue_class) {
            case 32:
            case 33:
            case 76:
            case 111:
            case 112:
            case 161:
            case 162:
            case 192:
            case 193:
            case 194:
            case 195:
            case 196:
            case 197:
            case 198:
            case 199:
                gll = 2.0;
                gtt = 1.0;
                break;
            case 80:
            case 81:
            case 82:
            case 83:
            case 84:
            case 85:
            case 86:
                // sinus node and surroundings
                gll = 0.5;
                gtt = 0.5;
                break;
            case 72:
                // crista terminalis
                gll = 7.7;
                gtt = 0.7;
                break;
            case 74:
            case 98:
            case 102:
            case 103:
                // pectinatae muscles, Bachman Bundle
                gll = 5.5;
                gtt = 2.75;
                break;
            case 104:
                gll = 1.1;
                gtt = 1.1;
                break;
            default:
                EXCEPTION("Unknown tissue class " << tissue_class << " at " << elementIndex);
        }

        mTensor(0,0) = gll;
        for (unsigned i = 1; i < DIM; i++)
            mTensor(i, i) = gtt;
    }
    
    c_matrix<double,DIM,DIM>& rCalculateModifiedConductivityTensor(unsigned elementIndex, const c_matrix<double,DIM,DIM>& rOriginalConductivity, unsigned domainIndex)
    {
        Element<DIM, DIM> *ele = pMesh->GetElement(elementIndex);
        if (ele->GetNumElementAttributes() > 0)
            ApplyTissueConductivity(elementIndex, (unsigned)ele->rGetElementAttributes()[0]);
        else
            mTensor.assign(rOriginalConductivity);

        if (conductivities.size() > 0)
            mTensor *= conductivities[elementIndex];

        return mTensor;
    }
};

template<unsigned DIM>
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
                CellMaleckar2008_baseFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                CellMaleckar2008_cAFFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                CellMaleckar2008_LA_1h2HzFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                CellMaleckar2008_RA_1h2HzFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                Cellcourtemanche_ramirez_nattel_1998_SRFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
                Cellcourtemanche_ramirez_nattel_1998_cAFFromCellMLCvodeOpt(noSolver, noStim).GetLookupTableCollection(),
        };

        for (int i = 0; i < (sizeof(tables) / sizeof(tables[0])); i++) {
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
        double stim_dur = GetDoubleOption("-stim_dur", 1.0);
        double stim_amp = GetDoubleOption("-stim_amp", 80000.0);
        LOG("\tduration : " << stim_dur << "ms");
        LOG("\tamplitude: " << stim_amp << "ms");

        LOG("\tcycles   : " << times.size());
        if (!times.empty()) {
            LOG("\tfirst    : " << times[0] << "ms");
            LOG("\tlast     : " << times.back() << "ms");
            if (user_defined) {
                std::stringstream ss;
                ss << times[0];
                for (unsigned i = 1; i < times.size(); i++)
                    ss << ", " << times[i];

                LOG("\ttimes    : " << ss.str());
            }
            else
                LOG("\tperiod   : " << period << "ms");
        }

        rStimTimes.insert( rStimTimes.end(), times.begin(), times.end() );
        return boost::shared_ptr<AbstractStimulusFunction>(new TimedStimulus(-stim_amp, stim_dur, times));
    }

    AtrialCellFactory<DIM> InitCellFactory(std::vector<double> &rStimTimes) {
        LOG("** CELLS **")
        std::string cellopt = "courtemanche_sr";
        if (CommandLineArguments::Instance()->OptionExists("-cell"))
            cellopt = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-cell");

        int cell_model;
        if (cellopt == "maleckar")
            cell_model = CellModel::MALECKAR;
        else if (cellopt == "maleckar_caf")
            cell_model = CellModel::MALECKAR_CAF;
        else if (cellopt == "maleckar_anna")
            cell_model = CellModel::MALECKAR_ANNA;
        else if (cellopt == "courtemanche_sr")
            cell_model = CellModel::COURTEMANCHE_SR;
        else if (cellopt == "courtemanche_caf")
            cell_model = CellModel::COURTEMANCHE_CAF;
        else
            EXCEPTION("Unknown Cell Model: " << cellopt);
        LOG("cell: " << cellopt);

        auto p_stim_sinus = InitStimulus("sinus", rStimTimes, 4, 0, 500);
        auto p_stim_extra = InitStimulus("extra", rStimTimes, 6, 400, 300);

        OverrideVoltageLookupRange();
        return AtrialCellFactory<DIM>(p_stim_sinus, p_stim_extra, cell_model);
    }

    void ApplyPerm(std::vector<unsigned int> &nodes, const std::vector<unsigned int> &permutation) {
        if (permutation.size() > 0) {
            for (unsigned i = 0; i < nodes.size(); i++)
                nodes[i] = permutation[nodes[i]];
        }

        std::sort(nodes.begin(), nodes.end());
    }

    void AddActivationMap(MonodomainProblem<DIM> *problem, const std::vector<double> &rStimTimes) {
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

    chaste::parameters::v2017_1::media_type GetFibreOrientation(std::string meshfile) {
        if (FileFinder(meshfile + ".ortho", RelativeTo::AbsoluteOrCwd).IsFile())
            return cp::media_type::Orthotropic;

        return cp::media_type::Axisymmetric;
    }

    MonodomainProblem<DIM> *InitProblem(AtrialCellFactory<DIM> *cell_factory, AtrialConductivityModifier<DIM> *conductivity_modifier)
    {
        CommandLineArguments* args = CommandLineArguments::Instance();
        HeartConfig* heartConfig = HeartConfig::Instance();
        MonodomainProblem<DIM>* problem;

        if (args->OptionExists("-svi"))
            heartConfig->SetUseStateVariableInterpolation(true);

        LOG("** PROBLEM **")
        if (args->OptionExists("-loaddir")) {
            std::string loaddir = args->GetStringCorrespondingToOption("-loaddir");
            LOG("loaddir: " << loaddir);
            problem = CardiacSimulationArchiver<MonodomainProblem<DIM> >::Load(loaddir);
        }
        else {
            std::string meshfile = args->GetStringCorrespondingToOption("-meshfile");
            LOG("meshfile: " << meshfile);
            heartConfig->SetMeshFileName(meshfile, GetFibreOrientation(meshfile));

            problem = new MonodomainProblem<DIM>(cell_factory);
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

        conductivity_modifier->SetMesh(&problem->rGetMesh());
        problem->GetTissue()->SetConductivityModifier(conductivity_modifier);

        return problem;
    }

    AtrialConductivityModifier<DIM> InitConductivities() {
		double base_cond = GetDoubleOption("-base_cond", 1.75);
		double anisotropy_ratio = GetDoubleOption("-ar", 9.21);
		LOG("\tbase conductivity: " << base_cond);
		LOG("\tanisotropy ratio : " << anisotropy_ratio);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(base_cond, base_cond/anisotropy_ratio, base_cond/anisotropy_ratio));

        std::vector<float> conductivities;
        if (CommandLineArguments::Instance()->OptionExists("-condmod")) {
            std::string path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-condmod");
            LOG("conductivities   : " << path);
            conductivities = ConductivityReader::ReadConductivities(FileFinder(path, RelativeTo::AbsoluteOrCwd));
        }

        return AtrialConductivityModifier<DIM>(conductivities);
    }

    void Save(MonodomainProblem<DIM> *problem) {
        if (CommandLineArguments::Instance()->OptionExists("-savedir"))
        {
            std::string savedir = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-savedir");
            LOG("savedir: " << savedir);
            CardiacSimulationArchiver<MonodomainProblem<DIM> >::Save(*problem, savedir);
        }
    }

    void WriteLog(OutputFileHandler out_dir)
    {
        if (!PetscTools::AmMaster())
            return;

        out_stream os = out_dir.OpenOutputFile("log.txt");
        *os << QutemuLog::GetLog();
        os->close();
    }

    void WritePermutation(OutputFileHandler out_dir, MonodomainProblem<DIM> *problem)
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
    void RunSimulation() throw(Exception)
    {
        double start_time = Timer::GetWallTime();
        SetSchemaLocations();

        COUT("Initializing");
        OutputFileHandler out_dir = InitOutput();
        InitTimesteps();
        std::vector<double> stim_times;
        AtrialCellFactory<DIM> cell_factory = InitCellFactory(stim_times);
        AtrialConductivityModifier<DIM> conductivity_modifier = InitConductivities();
        MonodomainProblem<DIM>* problem = InitProblem(&cell_factory, &conductivity_modifier);
        AddActivationMap(problem, stim_times);

        COUT("Solving");
        problem->Solve();
        Save(problem);

        HeartEventHandler::Headings();
        HeartEventHandler::Report();

        if (PetscTools::AmMaster())
            out_dir.FindFile("progress_status.txt").Remove();

        WritePermutation(out_dir, problem);

        LOG("finished: " << std::setprecision(3) << std::fixed << (Timer::GetWallTime() - start_time) << "s");
        WriteLog(out_dir);

        delete problem;
        COUT("Success");
    }
};

void GetNextLineFromStream(std::ifstream& rFileStream, std::string& rRawLine)
{
    bool line_is_blank;
    do
    {
        getline(rFileStream, rRawLine, '\n');
        if (rFileStream.eof())
            EXCEPTION("File contains incomplete data: unexpected end of file.");

        // Get rid of any comment
        rRawLine = rRawLine.substr(0, rRawLine.find('#',0));

        line_is_blank = (rRawLine.find_first_not_of(" \t",0) == std::string::npos);
    }
    while (line_is_blank);
}

bool Is2dMesh()
{
    if (!CommandLineArguments::Instance()->OptionExists("-meshfile"))
        return false;

    std::string nodeFileName = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-meshfile") + ".node";
    ifstream nodeFile(FileFinder(nodeFileName, RelativeTo::AbsoluteOrCwd).GetAbsolutePath().c_str());
    if (!nodeFile.is_open())
        return false;

    std::string buffer;
    GetNextLineFromStream(nodeFile, buffer);
    std::stringstream node_header_line(buffer);
    unsigned mNumNodes, dimension;
    node_header_line >> mNumNodes >> dimension;
    return dimension == 2;
}

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
        if (Is2dMesh())
            AtrialFibrosis<2>().RunSimulation();
        else
            AtrialFibrosis<3>().RunSimulation();
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
