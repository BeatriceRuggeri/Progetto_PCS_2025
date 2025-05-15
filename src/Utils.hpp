#pragma once

#include <iostream>
#include "PolyhedralMesh.hpp"

using namespace std;

namespace calcolo{
	double EdgeLength(const double& x1, const double& y1,const double& x2, const double& y2);
}

namespace PolyhedralLibrary
{

bool ImportMesh(PolyhedralMesh& mesh);

bool ImportCell0Ds(PolyhedralMesh& mesh);


bool ImportCell1Ds(PolyhedralMesh& mesh);

bool ImportCell2Ds(PolyhedralMesh& mesh);

bool CheckEdges(const PolyhedralMesh& mesh, const double& eps);

//bool CheckAreas(const PolygonalMesh& mesh, const double& eps);


}
