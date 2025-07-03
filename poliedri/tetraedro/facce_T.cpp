#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<int>> facce_tetraedro = {
			{0, 0, 2, 1, 1, 3, 0}, //A-B-C -- AC,BC,AB
			{1, 0, 1, 3, 0, 4, 2}, //A-B-D -- AB,BD,AD
			{2, 0, 3, 2, 2, 1, 5}, //A-C-D -- AD,AC,CD
			{3, 1, 2, 3, 5, 3, 4} //B-C-D -- CD,BC,BD
		};