#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<int>> spigoli_ottaedro = {
			{0, 0, 2}, // AC
			{1, 0, 3}, // AD
			{2, 0, 4}, // AE
			{3, 0, 5}, // AF
			{4, 1, 2}, // BC
			{5, 1, 3}, // BD
			{6, 1, 4}, // BE
			{7, 1, 5}, // BF
			{8, 2, 4}, // CE
			{9, 2, 5}, // CF
			{10, 3, 4}, // DE
			{11, 3, 5} // DF
		};