
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

//#include "../src/ICCCBDerivedCa.hpp"
#include "../src/imtiaz_2002d_noTstart_COR.hpp"
#include "../src/DummyCell.hpp"

using namespace std;

class ICCCellFactory : public AbstractCardiacCellFactory<3>
{

public:
    ICCCellFactory()
        : AbstractCardiacCellFactory<3>()
    {
    }

    //AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        double x = pNode->rGetLocation()[0];
        double y = pNode->rGetLocation()[1];
        double z = pNode->rGetLocation()[2];
        double r = 0.1;
        //ICCCBDerivedCa* cell = new ICCCBDerivedCa(mpSolver, mpZeroStimulus);
        Cellimtiaz_2002d_noTstart_CORFromCellML* cell = new Cellimtiaz_2002d_noTstart_CORFromCellML(mpSolver, mpZeroStimulus);
        cell->SetParameter("eta", 0.045);
        //cell->SetParameter("V_excitation", -60);
	      //cell->SetParameter("live_time", 10000);
        //cell->SetParameter("ode_time_step", 0.1);
        //cell->SetParameter("IP3Par", 0.00069);
        //cell->SetParameter("t_start", 600000);

        if (((x-0.0378)*(x-0.0378)+(y+1.336)*(y+1.336)+(z+2.305)*(z+2.305)) < r*r)
        {
            //cell->SetParameter("t_start", 0);
            cell->SetParameter("eta", 0.037);
        }
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
        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/rat_cm_32_32_8_lm_32_32_2.1", cp::media_type::Orthotropic);

        // Output visualization options, we ask for meshalyzer and cmgui
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        // Simulation settings
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(2000);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,250);

        // Output file/folder
        HeartConfig::Instance()->SetOutputDirectory("TestRatStomach3D_Monodomain_rat_cm_32_32_8_lm_32_32_2_2s_250ms_xi_imtiaz");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        ICCCellFactory cell_factory;
        MonodomainProblem<3> monodomain_problem(&cell_factory);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.000001, 0.003, 0.003));//(0.03, 3.4, 1.0));
        //HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.003, 0.003, 0.000001));//(0.03, 3.4, 1.0));

        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};
