#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<double>> vertici_ottaedro = {
			{0, 1.0, 0.0, 0.0}, // A
			{1, -1.0, 0.0, 0.0}, // B
			{2, 0.0, 1.0, 0.0}, // C
			{3, 0.0, -1.0, 0.0}, // D
			{4, 0.0, 0.0, 1.0}, // E
			{5, 0.0, 0.0, -1.0} // F
		};