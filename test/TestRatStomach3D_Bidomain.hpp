
#include <cxxtest/TestSuite.h>
#include "BidomainProblem.hpp"
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
        double r = 0.5;
	    CellICCBioPhy* cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
        cell->SetParameter("V_excitation", -60);
	      cell->SetParameter("live_time", 10000);
        cell->SetParameter("ode_time_step", 0.001);
        cell->SetParameter("IP3Par", 0.00069);
        cell->SetParameter("t_start", 600000);

        if (((x-0.0378)*(x-0.0378)+(y+1.336)*(y+1.336)+(z+2.305)*(z+2.305)) < r*r)
        {
          cell->SetParameter("t_start", 0);
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

        // Output file/folder
        HeartConfig::Instance()->SetOutputDirectory("TestRatStomach3D_Bidomain_test");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        // Output visualization options, we ask for meshalyzer and cmgui
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        // Simulation settings
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        //HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        //HeartConfig::Instance()->SetKSPPreconditioner("jacobi");
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(500);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.001,0.001,50);

        ICCCellFactory cell_factory;
        BidomainProblem<3> bidomain_problem(&cell_factory);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.003, 0.003, 0.000001));//(0.03, 3.4, 1.0));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.003, 0.003, 0.000001));//(0.03, 3.4, 1.0));

        bidomain_problem.Initialise();
        bidomain_problem.Solve();
    }
};
