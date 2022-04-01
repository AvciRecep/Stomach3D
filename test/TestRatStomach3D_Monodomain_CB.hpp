
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
private:
  	std::set<unsigned> setICCNode;
public:
    ICCCellFactory(std::set<unsigned> iccNodes)
        : AbstractCardiacCellFactory<3>(), setICCNode(iccNodes)

    {
        ReadLaplaceFile();
    }

    void ReadLaplaceFile()
    {
        std::ifstream inLaplaceInfo("projects/mesh/Stomach3D/rat_scaffold_16_16_2.1_laplace_longi_sw.txt");
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
        unsigned index = pNode->GetIndex();

        if (setICCNode.find(index) != setICCNode.end())
        {
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
              cell->SetParameter("V_excitation", -55);
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
        else
            return new CellDummyCellFromCellML(mpSolver, mpZeroStimulus);
    }
};

class TestRatStomach3D : public CxxTest::TestSuite
{
public:
    void TestStomach3D() //throw (Exception)
    {
        ///// Input file
        string fname = "rat_scaffold_16_16_2.1";
        HeartConfig::Instance()->SetMeshFileName("projects/mesh/Stomach3D/"+fname, cp::media_type::Orthotropic);

        ///// Simulation settings
        int sim_dur = 10000; // ms
        int write_freq = 1000; //ms
        HeartConfig::Instance()->SetSurfaceAreaToVolumeRatio(2000);
        HeartConfig::Instance()->SetUseAbsoluteTolerance(1e-3);
        HeartConfig::Instance()->SetCapacitance(3);
        HeartConfig::Instance()->SetSimulationDuration(sim_dur);  //ms.
        HeartConfig::Instance()->SetOdePdeAndPrintingTimeSteps(0.1,1,write_freq);

        ///// Output file/folder
        string out_path = "test_stomach3d_monodomain_"+fname+"_"+std::to_string(sim_dur)+"ms_"+std::to_string(write_freq)+"ms";
        //string out_add = "_"+std::to_string(c_fibre)+"_"+std::to_string(c_sheet)+"_"+std::to_string(c_normal)+"_ref10s_v1";
        string out_add = "_conductivity_test_v0";
        HeartConfig::Instance()->SetOutputDirectory(out_path+out_add);
        HeartConfig::Instance()->SetOutputFilenamePrefix("results");

        ///// Output visualization options
        HeartConfig::Instance()->SetVisualizeWithCmgui(false);
        HeartConfig::Instance()->SetVisualizeWithMeshalyzer(true);
        HeartConfig::Instance()->SetVisualizeWithVtk(true);

        ///// Cell factory
        TrianglesMeshReader<3,3> reader("projects/mesh/Stomach3D/"+fname);
        DistributedTetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        std::set<unsigned> iccNodes;
        // Iterating trough all Elements in the mesh and assigning attributes, conductivities and saving all ICC nodes
        for (DistributedTetrahedralMesh<3,3>::ElementIterator iter = mesh.GetElementIteratorBegin();
                        iter != mesh.GetElementIteratorEnd(); ++iter)
        {
            // Read Attributes
            double attribute = iter->GetAttribute();
            // Copy all nodes of the element to the elementIndexesICC list
            if (attribute == 1) // Check if ICC node
            {
                for(int j = 0; j<=3; ++j)
                {
                    iccNodes.insert(iter->GetNodeGlobalIndex(j));
                }
            }
        }
        ICCCellFactory cell_factory(iccNodes);

        ///// MonodomainProblem
        MonodomainProblem<3> monodomain_problem(&cell_factory);
        ///// conductivities
        double c_fibre = 0.01;
        double c_sheet = 0.10;
        double c_normal = 0.01;
        HeartConfig::Instance()->SetIntracellularConductivities(Create_c_vector(c_fibre, c_sheet, c_normal));
        ///// Solve
        monodomain_problem.Initialise();
        monodomain_problem.Solve();
    }
};
