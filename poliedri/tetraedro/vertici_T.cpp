#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<double>> vertici_tetraedro = {
			{0, 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}, // A
			{1, 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0)}, // B
			{2, -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0)}, // C
			{3, -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0)} // D
		};