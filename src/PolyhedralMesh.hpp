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
		
		vector<unsigned int> M0DMarker = {}; // marker vertici
		vector<unsigned int> M1DMarker = {}; // marker spigoli
		
		unsigned int p, q, b, c;
		unsigned int V; //num vertici
		unsigned int E; //num spigoli
		unsigned int nodo_i;
		unsigned int nodo_f;
		vector<unsigned int> cammino; //vettore che segna gli id del cammino 
	}
}