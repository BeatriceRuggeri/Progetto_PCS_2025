#pragma once

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary{
	struct PolyhedralMesh {
		vector<vector<double>> M0D;
		vector<vector<double>> M1D;
		vector<vector<double>> M2D;
		vector<double> M3D;
	}
}