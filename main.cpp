#include <iostream>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;


// cmake -S ./ -B ./Debug -DCMAKE_BUILD_TYPE="Debug"

int main()
{
	
	/*
	int Matrix PQ[5][2] = {{3, 3}, {3, 4}, {4, 3}, {5, 3}, {3, 5}};
	
	for (i= 0; i<=5; i++)
	{	
		if (PQ[i]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE1(P, Q, b, c)
		}	
		if (PQ[i][1]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE2 (P, Q, b, c)
		}
	}
	*/
	//trial
	
	PolyhedralLibrary::PolyhedralMesh mesh;
	
	bool success = ExportCell0Ds("Cell0Ds.txt",mesh.Vert_tetrahedron,3);
	
	if(!success){
		cerr<<"Export failed."<<endl;
		return 1;
	}
				
	
	
	
	return 0;
}
