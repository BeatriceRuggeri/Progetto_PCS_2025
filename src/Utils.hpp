#pragma once

#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;
/*
namespace calcolo{
	double EdgeLength(const double& x1, const double& y1,const double& x2, const double& y2);
}
*/
namespace PolyhedralLibrary
{

bool ExportMesh(PolyhedralMesh& mesh);

bool ExportCell0Ds(const std::string& outputFilePath, const double vertici[][3], int n);

//bool ExportCell0Ds(PolyhedralMesh& mesh);


//bool ExportCell1Ds(PolyhedralMesh& mesh);

//bool EXportCell2Ds(PolyhedralMesh& mesh);

//bool CheckEdges(const PolyhedralMesh& mesh, const double& eps);

//bool CheckAreas(const PolygonalMesh& mesh, const double& eps);


}
