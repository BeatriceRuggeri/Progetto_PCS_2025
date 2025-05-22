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
	
	bool success_edges_cube = ExportCell1Ds("C_Cell1Ds.txt",mesh.Cell1Ds);
	
	if(!success_edges_cube){
		cerr<<"Export cube's edges failed."<<endl;
		return 1;
	}
	
	/*
	bool success_tetrahedron = ExportCell0Ds("T_Cell0Ds.txt",mesh.Vert_tetrahedron,4);
	
	if(!success_tetrahedron){
		cerr<<"Export tetrahedron failed."<<endl;
		return 1;
	}
	
	
	bool success_cube = ExportCell0Ds("C_Cell0Ds.txt",mesh.Vert_cube,8);
	
	if(!success_cube){
		cerr<<"Export cube failed."<<endl;
		return 1;
	}
	
	bool success_octahedron = ExportCell0Ds("O_Cell0Ds.txt",mesh.Vert_octahedron,6);
	
	if(!success_octahedron){
		cerr<<"Export octahedron failed."<<endl;
		return 1;
	}
	
	bool success_dodecahedron = ExportCell0Ds("D_Cell0Ds.txt",mesh.Vert_dodecahedron,20);
	
	if(!success_dodecahedron){
		cerr<<"Export dodecahedron failed."<<endl;
		return 1;
	}
	
	bool success_icosahedron = ExportCell0Ds("I_Cell0Ds.txt",mesh.Vert_icosahedron,12);
	
	if(!success_icosahedron){
		cerr<<"Export icosahedron failed."<<endl;
		return 1;
	}
	*/			
	
	
	
	return 0;
}
