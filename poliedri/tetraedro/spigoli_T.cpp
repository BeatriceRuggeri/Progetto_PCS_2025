#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<int>> spigoli_tetraedro = {
			{0, 0, 1}, // AB
			{1, 0, 2}, // AC
			{2, 0, 3}, // AD
			{3, 1, 2}, // BC
			{4, 1, 3}, // BD
			{5, 2, 3} // CD
		};