#include <iostream>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main()
{
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
				
	
	return 0;
}
