
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
#include "AbstractConductivityModifier.hpp"
#include <OrthotropicConductivityTensors.hpp>
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

// *************************** CELL FACTORY ************************************* //
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
              if(info.x > x-0.001 && info.x < x+0.001  && info.y > y-0.001 && info.y < y+0.001 && info.z > z-0.001 && info.z < z + 0.001)
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
        cell->SetParameter("live_time", 12000);

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

// *************************** CONDUCTIVITY MODIFIER ************************************* //
class ConductivityModifier : public AbstractConductivityModifier<3,3>
{
    private:
        c_matrix<double, 3,3> mSpecialMatrix;
    public:
        ConductivityModifier()
            : AbstractConductivityModifier<3,3>(),
              mSpecialMatrix(zero_matrix<double>(3,3))
        {
        }

        c_matrix<double,3,3>& rCalculateModifiedConductivityTensor(unsigned elementIndex,
                                                                   const c_matrix<double,3,3>& rOriginalConductivity,
                                                                   unsigned domainIndex)
        {
            PRINT_4_VARIABLES(elementIndex,rOriginalConductivity(0,0),rOriginalConductivity(1,1),rOriginalConductivity(2,2));
            // Conductivities for ICC
            if (domainIndex == 0)
            {
                mSpecialMatrix(0,0) = 0.03;
                mSpecialMatrix(1,1) = 0.30;
                mSpecialMatrix(2,2) = 0.03;
            }
            return mSpecialMatrix;
        }
};

// *************************** SIMULATION ************************************* //
class TestRatStomach3D : public CxxTest::TestSuite
{

public:
    void TestStomach3D() //throw (Exception)
    {
        ///// Input file
        string fname = "rat_scaffold_nopylorus_noesophagus_16_16_1.1";
        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/"+fname, cp::media_type::Orthotropic);

        ///// Simulation settings
        int sim_dur = 1; // ms
        int write_freq = 1; //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(sim_dur);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,write_freq);
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.03, 0.30, 0.03));//(0.03, 3.4, 1.0));

        ///// Output file/folder
        string out_path = "test_stomach3d_monodomain_"+fname+"_"+std::to_string(sim_dur)+"ms_"+std::to_string(write_freq)+"ms";
        string out_add = "_conductivity_test_v2";
        HeartConfig::Instance()->SetOutputDirectory(out_path+out_add);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");
        HeartConfig::Instance()->SetOutputUsingOriginalNodeOrdering(true);

        ///// Output visualization options
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        ///// Cell factory
        ICCCellFactory cell_factory;

        ///// MonodomainProblem
        MonodomainProblem<3> monodomain_problem(&cell_factory);
        monodomain_problem.Initialise();

        ///// Set conductivities
        //MonodomainTissue<3>* p_monodomain_tissue = monodomain_problem.GetMonodomainTissue();
        //ConductivityModifier modifier; // Initialise Conductivity Modifier
        //p_monodomain_tissue->SetConductivityModifier(&modifier);

        ///// Solve
        monodomain_problem.Solve();
    }
};
