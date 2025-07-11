#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<double>> vertici_icosaedro = {
			{0, 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, // A
			{1, 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, // B
			{2, 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))}, // C
			{3, 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))}, // D
			{4, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0}, // E
			{5, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0}, // F
			{6, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0}, // G
			{7, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0}, // H
			{8, (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, // I
			{9, (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))}, //L
			{10, -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, // M
			{11, -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))} //N
		};