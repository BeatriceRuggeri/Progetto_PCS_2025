#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void triangolazione_I(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshOutput)
	{
		PolyhedralMesh meshTriangolazione;

		vector<int> dimension = topological_Calc_I(q, b, c);
		vector<int> dimensionDuplicated = duplicates_Calc(q, b, c, dimension);
		tri_build_I(mesh, meshTriangolazione, b, c, dimensionDuplicated);
		v_duplicates_rm(meshTriangolazione);
    	l_duplicates_rm(meshTriangolazione);
    	meshClean(meshTriangolazione, meshOutput, dimension);
    	Assemble_3D(meshOutput);
    }
    
    void triangolazione_dual(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshDual)
	{
		PolyhedralMesh meshTriangolazione;
		PolyhedralMesh meshOutput;
		
		vector<int> dimension = topological_Calc_I(q, b, c);
		vector<int> dimensionDuplicated = duplicates_Calc(q, b, c, dimension);
		tri_build_I(mesh, meshTriangolazione, b, c, dimensionDuplicated);
		v_duplicates_rm(meshTriangolazione);
    	l_duplicates_rm(meshTriangolazione);
    	meshClean(meshTriangolazione, meshOutput, dimension);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildMap_LF(meshOutput);
    	dual_calc(meshOutput, meshDual, edgeToFacesMap);
		Assemble_3D(meshDual);
    }
    
    void triangolazione_II(const int q, const int b, PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolazione_II)
	{
		PolyhedralMesh meshTriangolazione;
		PolyhedralMesh meshOutput;
		
		vector<int> dimension = topological_Calc_I(q, b, 0);
		vector<int> dimensionDuplicated = duplicates_Calc(q, b, 0, dimension);
		tri_build_I(mesh, meshTriangolazione, b, 0,  dimensionDuplicated);
		v_duplicates_rm(meshTriangolazione);
    	l_duplicates_rm(meshTriangolazione);
    	meshClean(meshTriangolazione, meshOutput, dimension);
    	
    	vector<int> dimension2 = topological_Calc_II(b, q);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildMap_LF(meshOutput);
		tri_build_II(meshOutput, meshTriangolazione_II, dimension2, edgeToFacesMap);
		Assemble_3D(meshTriangolazione_II);
    }

	void Assemble_3D(PolyhedralMesh& meshTriangolazione){
		
		meshTriangolazione.Cell3DsId = {0};
		meshTriangolazione.NumCells0Ds = meshTriangolazione.Cell0DsId.size();
		meshTriangolazione.NumCells1Ds = meshTriangolazione.Cell1DsId.size();
		meshTriangolazione.NumCells2Ds = meshTriangolazione.Cell2DsId.size();

		for (unsigned int i = 0; i < meshTriangolazione.Cell0DsId.size(); i++){
			meshTriangolazione.Cell3DsVertices.push_back(meshTriangolazione.Cell0DsId[i]);
		}
		for (unsigned int i = 0; i < meshTriangolazione.Cell1DsId.size(); i++){
			meshTriangolazione.Cell3DsEdges.push_back(meshTriangolazione.Cell1DsId[i]);

		}
		for (unsigned int i = 0; i < meshTriangolazione.Cell2DsId.size(); i++){
			meshTriangolazione.Cell3DsFaces.push_back(i);
		}
	}
	

}