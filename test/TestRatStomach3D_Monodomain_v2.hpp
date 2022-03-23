
#include <cxxtest/TestSuite.h>
#include "MonodomainProblem.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractStimulusFactory.hpp"
#include "ElectrodesStimulusFactory.hpp"

#include "ChasteNodesList.hpp"
#include "ChasteCuboid.hpp"
#include "AbstractChasteRegion.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"
#include "ArchiveOpener.hpp"
#include "AbstractCardiacProblem.hpp"

#include "../src/CellICCBioPhy.hpp"
#include "../src/DummyCell.hpp"

using namespace std;

struct coordinateV_st
{
    double x;
    double y;
    double z;
    double V;
};

class ICCCellFactory : public AbstractCardiacCellFactory<3>
{
private:
    std::vector<coordinateV_st> LaplaceInfo;

public:
    ICCCellFactory()
        : AbstractCardiacCellFactory<3>()
    {
        ReadLaplaceFile();
    }

    void ReadLaplaceFile()
    {
        std::ifstream inLaplaceInfo("projects/mesh/Stomach3D/rat_scaffold_nopylorus_noesophagus_16_16_1.1_laplace_longi_sw.txt");
        if(!inLaplaceInfo)
        {
          EXCEPTION("Reading laplace solution error");
        }
        std::string line;
        coordinateV_st lapInfo;

        while(std::getline(inLaplaceInfo, line))
        {
          stringstream cordinateLap(line);
          cordinateLap >> lapInfo.x >> lapInfo.y >> lapInfo.z >> lapInfo.V;
          LaplaceInfo.push_back(lapInfo);
        }
    }

    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];
        coordinateV_st info;
        int counter = 0;
        double V_val = 0;
        for(std::vector<coordinateV_st>::iterator itr = LaplaceInfo.begin(); itr!=LaplaceInfo.end();itr++)
        {
            info = *itr;
            if(info.x > x-0.001 && info.x < x+0.001  && info.y > y-0.001 && info.y < y+0.001 && info.z >
                   z-0.001 && info.z < z + 0.001)
            {
                counter++;
                V_val = info.V;
                break;
            }
        }
        if (counter != 1)
        {
          PRINT_4_VARIABLES(x,y,z, counter);
          EXCEPTION("Coordinates not found in Laplace file");
        }


	      CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
        cell->SetParameter("V_excitation", -60);
	      cell->SetParameter("live_time", 8000);

        if (V_val > 97)
        {
            return new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
        }
        double r = 0.5;
        if (((x-0)*(x-0)+(y+1.46)*(y+1.46)+(z+2.7)*(z+2.7)) < r*r)
            cell->SetParameter("t_start", 0);
        else
            cell->SetParameter("t_start", 600000);

        cell->SetParameter("ode_time_step", 0.1);
        cell->SetParameter("IP3Par", 0.00069);
        return cell;
    }
};

class TestRatStomach3D : public CxxTest::TestSuite
{
public:
    void TestStomach3D() //throw (Exception)
    {
        // Input file
        //HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/rat_cm_lm_16_16_1", cp::media_type::Orthotropic);
        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/rat_scaffold_nopylorus_noesophagus_16_16_1.1", cp::media_type::Orthotropic);

        // Output visualization options, we ask for meshalyzer and cmgui
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        // Simulation settings
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(20000);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,250);

        // Output file/folder
        HeartConfig::Instance()->SetOutputDirectory("test_stomach3d_monodomain_20s_250ms_0.03_0.30_0.03_rat_scaffold_nopylorus_noesophagus_16_16_1.1");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        ICCCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        //HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.000005, 0.075, 0.0005));//(0.03, 3.4, 1.0));
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.03, 0.30, 0.03));

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};
