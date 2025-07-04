#pragma once


#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {
	
	struct PolyhedralMesh{
		
		vector<unsigned int> Cell0DsId = {};
		Eigen::MatrixXd Cell0DsCoordinates = {};
		vector<vector<unsigned int>> Cell0DsFlag= {};
		vector<unsigned int> Cell0DsMarker = {};
		
		vector<unsigned int> Cell1DsId = {}; 
		MatrixXi Cell1DsExtrema = {}; 
		
		
		vector<unsigned int> Cell1DsFlag= {}; 
		vector<bool> Cell1DsOriginalFlag = {}; // !!!!!!!
		vector<unsigned int> Cell1DsMarker = {};    
		
		vector<unsigned int> Cell2DsId = {}; 
		vector<vector<unsigned int>> Cell2DsVertices = {}; 
		vector<vector<unsigned int>> Cell2DsEdges = {}; 
		
		unsigned int Cell3DsId = 0; 
		unsigned int NumCells0Ds = 0;
		unsigned int NumCells1Ds = 0; 
		unsigned int NumCells2Ds = 0; 
		vector<unsigned int> Cell3DsVertices = {}; 
		vector<unsigned int> Cell3DsEdges = {}; 
		vector<unsigned int> Cell3DsFaces = {}; 

   };

}
