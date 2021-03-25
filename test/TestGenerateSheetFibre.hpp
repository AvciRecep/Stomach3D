
#ifndef TestGenerateSheetFibre_HPP_
#define TestGenerateSheetFibre_HPP_

#include <cxxtest/TestSuite.h>

#include "UblasIncludes.hpp"
/* This class represents the mesh internally. */
#include "TetrahedralMesh.hpp"
/* These are used to specify boundary conditions for the PDEs. */
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"
/* This class helps us deal with output files. */
#include "OutputFileHandler.hpp"
/* The following header must be included in every test that uses PETSc. Note that it
 * cannot be included in the source code. */
#include "PetscSetupAndFinalize.hpp"

// My includes
#include <math.h>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include <Eigen/Geometry>

#define PI 3.14159265

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::AngleAxisd;
using Eigen::Vector3d;
using Eigen::Quaterniond;
#include "VtkMeshWriter.hpp"

class TestFiberSheetHR : public CxxTest::TestSuite
{

private:

public:

    void TestCalculateSheetHR()
    {
        TrianglesMeshReader<3,3> mesh_reader("/media/hpc/CHASTE2017/chaste-source/projects/FibresHSB016/src/hsb16_tesel_8.1");
        /* Now declare a tetrahedral mesh with the same dimensions... */
        TetrahedralMesh<3,3> mesh;
        /* ... and construct the mesh using the mesh reader. */
        mesh.ConstructFromMeshReader(mesh_reader);

        std::ifstream gradijk("/media/hpc/CHASTE/2017-simulation-results/TestLaplaceHSB016/grad_longi.txt");
        if (!gradijk)
        {
            cout << "There was a problem opening laplace gradient for reading " << endl;
        }

        std::ifstream gradNormalijk("/media/hpc/CHASTE/2017-simulation-results/TestLaplaceHSB016/grad_circum.txt");
        if (!gradNormalijk)
        {
            cout << "There was a problem opening laplace gradient normal for reading " << endl;
        }
        
        OutputFileHandler output_file_handler("HSB016Full");
        out_stream p_file = output_file_handler.OpenOutputFile("stomach.ortho");
        double x, y, z;
        Vector3d fibre;
        Vector3d sheet;
        Vector3d normal;
        string line;
        std::vector<c_vector<double, 3u> > fibre_directions;
        std::vector<c_vector<double, 3u> > sheet_directions;
        std::vector<c_vector<double, 3u> > cross_directions;	
        while (std::getline(gradijk, line))
        {
            std::stringstream lapStream(line);
            lapStream >> x >> y >> z;
            fibre(0) = x;
            fibre(1) = y;
            fibre(2) = z;
            fibre = fibre.normalized();
            if (((fibre.array() != fibre.array())).all())
            {
                fibre(0) = 1;
                fibre(1) = 0;
                fibre(2) = 0;
            }
            
            c_vector<double, 3u> fibre_direction;
            fibre_direction[0] = fibre(0);
            fibre_direction[1] = fibre(1);
            fibre_direction[2] = fibre(2);
            fibre_directions.push_back(fibre_direction);
            /*********************************************/
            /* normal directions */
            /*********************************************/
            std::getline(gradNormalijk, line);
            std::stringstream lapNormStream(line);
            lapNormStream >> x >> y >> z;
            normal(0) = x;
            normal(1) = y;
            normal(2) = z;
            normal = normal.normalized();
            if (((normal.array() != normal.array())).all())
            {
                normal(0) = 0;
                normal(1) = 0;
                normal(2) = 1;
            }
            double dot_product = fibre.dot(normal);
            Vector3d temp;
            if( fabs(dot_product) != 0)
            {
                /*********************************************/
                /* Fix orthoganality */
                /*********************************************/
                temp = fibre.cross(normal);
                normal = temp.cross(fibre).normalized();
            }
            c_vector<double, 3u> normal_direction;
            normal_direction[0] = normal(0);
            normal_direction[1] = normal(1);
            normal_direction[2] = normal(2);
            cross_directions.push_back(normal_direction);
            
            /*********************************************/
            /* sheet directions */
            /*********************************************/
            sheet = fibre.cross(normal);
            c_vector<double, 3u> sheet_direction;
            sheet_direction[0] = sheet(0);
            sheet_direction[1] = sheet(1);
            sheet_direction[2] = sheet(2);
            sheet_directions.push_back(sheet_direction);
            /*********************************************/
            /* Write CHASTE ortho file */
            /*********************************************/
            if (((fibre.array() != fibre.array())).all() ||((normal.array() != normal.array())).all() ||((sheet.array() != sheet.array())).all() )
            {
                fibre(0) = 1;
                fibre(1) = 0;
                fibre(2) = 0;
                sheet(0) = 0;
                sheet(1) = 1;
                sheet(2) = 0;
                normal(0) = 0;
                normal(1) = 0;
                normal(2) = 1;
            }
            (*p_file) << fibre(0) << " " << fibre(1) << " " << fibre(2) << " "
                      << sheet(0) << " " << sheet(1) << " " << sheet(2) << " "
                      << normal(0) << " " << normal(1) << " " << normal(2) << "\n";
        }
        p_file->close();
        VtkMeshWriter<3u, 3u> mesh_writer("HSB016Full", "mesh", false);
        mesh_writer.AddCellData("Fibre Direction", fibre_directions);
        mesh_writer.AddCellData("Sheet Direction", sheet_directions);
        mesh_writer.AddCellData("Cross Direction", cross_directions);
        mesh_writer.WriteFilesUsingMesh(mesh);
    }
};

#endif 
