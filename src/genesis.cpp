#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary{
	
void tetraedro_gen(PolyhedralMesh& mesh) {
		
		double sr_3 =sqrt(3.0)/3.0;
		
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 4);
	
		mesh.Cell0DsCoordinates(0,0) = sr_3;   mesh.Cell0DsCoordinates(1,0) = sr_3;   mesh.Cell0DsCoordinates(2,0) = sr_3;  
		mesh.Cell0DsCoordinates(0,1) = -sr_3;  mesh.Cell0DsCoordinates(1,1) = -sr_3;  mesh.Cell0DsCoordinates(2,1) = sr_3; 
		mesh.Cell0DsCoordinates(0,2) = -sr_3;  mesh.Cell0DsCoordinates(1,2) = sr_3;   mesh.Cell0DsCoordinates(2,2) = -sr_3;
		mesh.Cell0DsCoordinates(0,3) = sr_3;   mesh.Cell0DsCoordinates(1,3) = -sr_3;  mesh.Cell0DsCoordinates(2,3) = -sr_3;
		
		mesh.Cell0DsId = {0,1,2,3};
		
//*********************************************************************************************************
		mesh.Cell1DsId = {0,1,2,3,4,5};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(6, 2);
		mesh.Cell1DsExtrema <<
			0, 1,
			1, 2, 
			2, 0,
			0, 3, 
			3, 1, 
			2, 3;
			
//*********************************************************************************************************
		mesh.Cell2DsId = {0,1,2,3};
	
		mesh.Cell2DsVertices = {
			{0,1,2},
			{0,3,1},
			{1,3,2},
			{2,3,0}
		};
	
		mesh.Cell2DsEdges = {
			{0,1,2},
			{3,4,0},
			{4,5,1},
			{5,3,2}
		};
	
//*********************************************************************************************************
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0,1,2,3};
		mesh.Cell3DsEdges = {0,1,2,3,4,5};
		mesh.Cell3DsFaces = {0,1,2,3};
	
	}
	
//*********************************************************************************************************

void icosaedro_gen(PolyhedralMesh& mesh) {
		
		double rad = 1.0;
		double phi = (1.0 + sqrt(5.0))/2.0;
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 12);
		
		double norm = sqrt(10.0+2.0*sqrt(5.0))/2.0;
		mesh.Cell0DsCoordinates(0, 0) = 0.0; mesh.Cell0DsCoordinates(1, 0) = (rad)/norm; mesh.Cell0DsCoordinates(2, 0) = phi/norm;
		mesh.Cell0DsCoordinates(0, 1) = 0.0; mesh.Cell0DsCoordinates(1, 1) = -rad/norm; mesh.Cell0DsCoordinates(2, 1) = phi/norm;
		mesh.Cell0DsCoordinates(0, 2) = 0.0; mesh.Cell0DsCoordinates(1, 2) = (rad)/norm; mesh.Cell0DsCoordinates(2, 2) = -phi/norm;
		mesh.Cell0DsCoordinates(0, 3) = 0.0; mesh.Cell0DsCoordinates(1, 3) = -(rad)/norm; mesh.Cell0DsCoordinates(2, 3) = -phi/norm;
		mesh.Cell0DsCoordinates(0, 4) = (rad)/norm; mesh.Cell0DsCoordinates(1, 4) = phi/norm; mesh.Cell0DsCoordinates(2, 4) = 0.0;
		mesh.Cell0DsCoordinates(0, 5) = -(rad)/norm; mesh.Cell0DsCoordinates(1, 5) = phi/norm; mesh.Cell0DsCoordinates(2, 5) = 0.0;
		mesh.Cell0DsCoordinates(0, 6) = (rad)/norm; mesh.Cell0DsCoordinates(1, 6) = -phi/norm; mesh.Cell0DsCoordinates(2, 6) = 0.0;
		mesh.Cell0DsCoordinates(0, 7) = -(rad)/norm; mesh.Cell0DsCoordinates(1, 7) = -phi/norm; mesh.Cell0DsCoordinates(2, 7) = 0.0;
		mesh.Cell0DsCoordinates(0, 8) = phi/norm; mesh.Cell0DsCoordinates(1, 8) = 0.0; mesh.Cell0DsCoordinates(2, 8) = (rad)/norm;
		mesh.Cell0DsCoordinates(0, 9) = -phi/norm; mesh.Cell0DsCoordinates(1, 9) = 0.0; mesh.Cell0DsCoordinates(2, 9) = (rad)/norm;
		mesh.Cell0DsCoordinates(0, 10) = phi/norm; mesh.Cell0DsCoordinates(1, 10) = 0.0; mesh.Cell0DsCoordinates(2, 10) = -(rad)/norm;
		mesh.Cell0DsCoordinates(0, 11) = -phi/norm; mesh.Cell0DsCoordinates(1, 11) = 0.0; mesh.Cell0DsCoordinates(2, 11) = -(rad)/norm;
	
		mesh.Cell0DsId = {0,1,2,3,4,5,6,7,8,9,10,11};
//*********************************************************************************************************
	
		mesh.Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
					 20,21,22,23,24,25,26,27,28,29};
					 
		
		mesh.Cell1DsExtrema = MatrixXi::Zero(30, 2);
		mesh.Cell1DsExtrema << //ends, LASCIATI COMMENTI
		  0, 1,   // 0
		  0, 4,   // 1
		  0, 5,   // 2
		  0, 8,   // 3
		  0, 9,   // 4
		  1, 6,   // 5
		  1, 7,   // 6
		  1, 8,   // 7
		  1, 9,   // 8
		  2, 3,   // 9
		  2, 10,   //10
		  2, 4,   //11
		  3, 6,   //12
		  3, 7,   //13
		  3, 10,  //14
		  3, 11,  //15
		  4, 5,   //16
		  9, 7,   //17
		  4, 8,   //18
		  4, 10,  //19
		  5, 2,   //20
		  5, 9,   //21
		  5, 11,  //22
		  6, 8,   //23
		  6, 7,   //24
		  6, 10,  //25
		  7, 11,  //26
		  8, 10,  //27
		  9, 11,  //28
		  2, 11;  //29

//*********************************************************************************************************
//faces
		mesh.Cell2DsId = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
		
		mesh.Cell2DsVertices = {
			{0,1,9},
			{0,9,5},
			{0,5,4},
			{0,4,8},
			{0,8,1},  
			{1,6,8}, 
			{8,6,10}, 
			{10,8,4}, 
			{10,4,2}, 
			{4,2,5}, 
			{5,2,11}, 
			{5,11,9}, 
			{9,7,11}, 
			{1,9,7}, 
			{1,7,6}, 
			{7,11,3}, 
			{11,3,2}, 
			{2,3,10}, 
			{10,6,3},
			{3,6,7}
		};
//*********************************************************************************************************

		
		mesh.Cell2DsEdges = {
			{0,8,4},     
			{4,21,2},    
			{2,16,1},    
			{1,18,3},    
			{3,7,0},     
			{5,23,7},    
			{23,25,27},  
			{27,18,19}, 
			{19,11,10}, 
			{11,20,16},  
			{20,29,22},  
			{22,28,21},  
			{17,26,28},  
			{8,17,6},    
			{6,24,5},    
			{26,15,13},   
			{15,9,29},  
			{9,14,10},  
			{25,12,14},  
			{12,24,13}    

		};
		//
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0,1,2,3,4,5,6,7,8,9,10,11};
		mesh.Cell3DsEdges = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
		mesh.Cell3DsFaces = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
		
	}
//*********************************************************************************************************

				
void octaedro_gen(PolyhedralMesh& mesh) {
		
		double oct_rad = 1.0; 
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 6);
	
		mesh.Cell0DsCoordinates(0,0) = oct_rad;   mesh.Cell0DsCoordinates(1,0) = 0.0;  mesh.Cell0DsCoordinates(2,0) = 0.0;  
		mesh.Cell0DsCoordinates(0,1) = -oct_rad;  mesh.Cell0DsCoordinates(1,1) = 0.0;  mesh.Cell0DsCoordinates(2,1) = 0.0;  
		mesh.Cell0DsCoordinates(0,2) = 0.0;  mesh.Cell0DsCoordinates(1,2) = oct_rad;   mesh.Cell0DsCoordinates(2,2) = 0.0;  
		mesh.Cell0DsCoordinates(0,3) = 0.0;  mesh.Cell0DsCoordinates(1,3) = -oct_rad;  mesh.Cell0DsCoordinates(2,3) = 0.0;  
		mesh.Cell0DsCoordinates(0,4) = 0.0;  mesh.Cell0DsCoordinates(1,4) = 0.0;  mesh.Cell0DsCoordinates(2,4) = oct_rad;   
		mesh.Cell0DsCoordinates(0,5) = 0.0;  mesh.Cell0DsCoordinates(1,5) = 0.0;  mesh.Cell0DsCoordinates(2,5) = -oct_rad;  
	
		mesh.Cell0DsId = {0,1,2,3,4,5};
		
//*********************************************************************************************************
//edges
		mesh.Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(12,2);
		mesh.Cell1DsExtrema <<
			0, 2,  
			2, 1, 
			1, 3,  
			3, 0,  
			0, 4,  
			2, 4,  
			4, 1,  
			4, 3,  
			1, 5,  
			3, 5,  
			0, 5,  
			2, 5;  
			
//*********************************************************************************************************
//faces
		mesh.Cell2DsId = {0,1,2,3,4,5,6,7};

		mesh.Cell2DsVertices = {
			{0,2,4}, 
			{0,4,3}, 
			{3,4,1},  
			{1,2,4}, 
			{2,5,0}, 
			{2,1,5},  
			{1,5,3},  
			{0,3,5}  
		};
	
		mesh.Cell2DsEdges = {
			{0,5,4},  
			{4,7,3},  
			{7,6,2}, 
			{1,5,6},  
			{11,10,0}, 
			{1,8,11}, 
			{8,9,2},  
			{3,9,10}  
		};
	
//polyhedral: all together
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0,1,2,3,4,5};
		mesh.Cell3DsEdges = {0,1,2,3,4,5,6,7,8,9,10,11};
		mesh.Cell3DsFaces = {0,1,2,3,4,5,6,7};
		
	}

}