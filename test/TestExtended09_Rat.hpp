#ifndef TESTEXTENDED09_HPP_
#define TESTEXTENDED09_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "../src/CellICCBioPhy.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "AbstractStimulusFactory.hpp"
#include "ElectrodesStimulusFactory.hpp"
#include "ExtendedBidomainProblem.hpp"
#include "ChasteNodesList.hpp"
#include "ChasteCuboid.hpp"
#include "AbstractChasteRegion.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "CorriasBuistSMCModified.hpp"
#include "../src/DummyCell.hpp"
#include "Debug.hpp"
#include "ArchiveOpener.hpp"
#include "AbstractCardiacProblem.hpp"


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
        std::ifstream inLaplaceInfo("projects/Stomach3D/src/rat_16_16_1_linear_sol_longi.txt");
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

    //AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {

        //double x = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[0];
      	//double y = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[1];
        //double z = this->GetMesh()->GetNode(nodeIndex)->rGetLocation()[2];
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
      	CellICCBioPhy* cell = 0;
      	cell = new CellICCBioPhy(mpSolver, mpZeroStimulus);
      	cell->SetParameter("V_excitation", -55);
      	cell->SetParameter("live_time", 12000);

        if (V_val > 50)
        {
            return new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
        }

        //ChastePoint<3> centre(7.24767, -2.34362, -1.79832); // for human stomach H09ext file
        ChastePoint<3> centre(0.03399, -1.038, -1.824); // for rat stomach mesh.
        ChastePoint<3> radii (0.1, 0.1, 0.1);
        ChasteEllipsoid<3> ellipseRegion(centre, radii);
        ChastePoint<3> myPoint(x, y, z);

        if(ellipseRegion.DoesContain(myPoint))
            cell->SetParameter("t_start", 0);
        else
            cell->SetParameter("t_start",500000);

        cell->SetParameter("ode_time_step",0.1);
        cell->SetParameter("IP3Par", 0.0006);
        return cell;
    }
};

class SMCCellFactory : public AbstractCardiacCellFactory<3>
{
public:
    SMCCellFactory() : AbstractCardiacCellFactory<3>()
    {
    }

    //AbstractCardiacCell* CreateCardiacCellForTissueNode(unsigned nodeIndex)
    AbstractCardiacCell* CreateCardiacCellForTissueNode(Node<3>* pNode)
    {
        CorriasBuistSMCModified *cell;
        cell = new CorriasBuistSMCModified(mpSolver, mpZeroStimulus);

        cell->SetFakeIccStimulusPresent(false);//it will get it from the real ICC, via gap junction
        return cell;
    }
};


class TestExtended60Stomach : public CxxTest::TestSuite
{
public:

    void TestExtendedStomach() //throw (Exception)
    {
        Extended();
        //ExtendedRemaining();
    }

private:
    void ExtendedRemaing() //throw (Exception)
    {
        FileFinder archive_dir("Extended", RelativeTo::ChasteTestOutput);
        std::string archive_file = "extended.arch";
        ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();
        AbstractCardiacProblem<3,3,3> *p_problem;
        (*p_arch) >> p_problem;
        TS_ASSERT_EQUALS(HeartConfig::Instance()->GetSimulationDuration(),5000);
        ExtendedBidomainProblem<3>* p_extended_problem = dynamic_cast<ExtendedBidomainProblem<3>*>(p_problem);
        unsigned size_first_cell, size_second_cell, size_extrastim, global_size_first_cell, global_size_second_cell, global_size_extrastim;
        for (unsigned proc = 0; proc < PetscTools::GetNumProcs(); proc++)
        {
            size_first_cell = p_problem->GetTissue()->rGetCellsDistributed().size();
            size_second_cell = p_extended_problem->GetExtendedBidomainTissue()->rGetSecondCellsDistributed().size();
            size_extrastim = p_extended_problem->GetExtendedBidomainTissue()->rGetExtracellularStimulusDistributed().size();
        }
        int mpi_ret_1 = MPI_Allreduce(&size_first_cell, &global_size_first_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        int mpi_ret_2 = MPI_Allreduce(&size_second_cell, &global_size_second_cell, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);
        int mpi_ret_3 = MPI_Allreduce(&size_extrastim, &global_size_extrastim, 1, MPI_UNSIGNED, MPI_SUM, PETSC_COMM_WORLD);

        assert(mpi_ret_1 == MPI_SUCCESS);
        assert(mpi_ret_2 == MPI_SUCCESS);
        assert(mpi_ret_3 == MPI_SUCCESS);

        HeartConfig::Instance()->SetSimulationDuration(20000);
        p_problem->Solve();
        delete p_problem;
    }

private:
    void Extended() //throw (Exception)
    {
        HeartConfig::Instance()->Reset();
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,1000);
        HeartConfig::Instance()->SetKSPSolver("gmres");
        HeartConfig::Instance()->SetUseAbsoluteTolerance(2e-1);
        HeartConfig::Instance()->SetKSPPreconditioner("jacobi");

        HeartConfig::Instance()->SetSimulationDuration(20000);  //ms.

        // Output visualization options, we ask for meshalyzer and cmgui
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(false);

        HeartConfig::Instance()->SetOutputDirectory("Stomach3D_rat_16_16_1_dt1s_20s_000003_0001_5");
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/rat_16_16_1.1", cp::media_type::Orthotropic);

        ICCCellFactory tissueICCInfo;
        SMCCellFactory tissueSMCInfo;

        ExtendedBidomainProblem<3> extended_problem(&tissueICCInfo, &tissueSMCInfo);

        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(0.000003,0.0001,0.5));//(0.03, 3.4, 1.0));
        HeartConfig::Instance()->SetExtracellularConductivities(Create_c_vector(0.000003,0.0001,0.5));//(0.03, 3.4, 1.0));

        double Am_icc = 2000.0;
        double Am_smc = 1000.0;
        double Am_gap = 1.0;
        double Cm_icc = 2.5;
        double Cm_smc = 1.0;
        double G_gap = 20.0;

        extended_problem.SetExtendedBidomainParameters(Am_icc,Am_smc, Am_gap, Cm_icc, Cm_smc, G_gap);
        extended_problem.SetIntracellularConductivitiesForSecondCell(Create_c_vector(0.000003,0.0001,0.5));//(0.02, 3.4, 1.0));

        extended_problem.Initialise();
        extended_problem.Solve();

        // CheckPoint the simulation
        FileFinder archive_dir("Extended", RelativeTo::ChasteTestOutput);
        std::string archive_file = "extended.arch";
        ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
        boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

        AbstractCardiacProblem<3,3,3>* const p_extended_problem = &extended_problem;
        (*p_arch) & p_extended_problem;

        HeartEventHandler::Headings();
        HeartEventHandler::Report();
    }
};
#endif
