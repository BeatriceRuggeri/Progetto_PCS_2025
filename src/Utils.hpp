
#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>

using namespace PolyhedralLibrary;
using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary{
//tanto questo sono l'impostazioni
	
void swap_val(int& p, int& q);
vector<int> topological_Calc_I(const int q, const int b, const int c);
vector<int> duplicates_Calc(const int q, const int b, const int c, const vector<int>& dim);
vector<int> topological_Calc_II(const int b, const int q);
	
//cleaning 
void v_duplicates_rm(PolyhedralMesh& meshTriangulated);
void l_duplicates_rm(PolyhedralMesh& meshTriangulated);
void meshClean(const PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension);	
	
	
	
//gen dati:
void tetraedro_gen(PolyhedralMesh& mesh);
void octaedro_gen(PolyhedralMesh& mesh);
void icosaedro_gen(PolyhedralMesh& mesh);
	
//Quali input?-> cambiamo a se p==3
	
void triangolazione_I(const int q,const int b,const int c,PolyhedralMesh& mesh,PolyhedralMesh& meshOutput);
//****************************** utili
void edge_correction(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangolazione, unsigned int& id_edge, const unsigned int id_triangle); //correct those missing ones here
unsigned int vertex_correction(const Vector3d& coord, PolyhedralMesh& meshTriangolazione, unsigned int& i1);
unsigned int edge_correction_II(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& i2);
Vector3d baricenter_F(const PolyhedralMesh& meshTriangolazione, const unsigned int edgeId, const unsigned int currentFaceId, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);
void tri_build_I(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolazione, const unsigned int b, const unsigned int c, const vector<int>& dimension);
vector<unsigned int> normalization_L(const vector<unsigned int>& face_edges);
vector<unsigned int> norm_cyclic(const vector<unsigned int>& current_edges);
   
//********************************************************

void triangolazione_II(const int q, const int b, PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolazione_II);
void face_F(const vector<unsigned int>& new_face_vertices, const vector<unsigned int>& new_face_edges, PolyhedralMesh& meshTriangulated, unsigned int& k3);
void tri_build_II(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolazione, const vector<int>& dimension, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edge_faces_Map);
void triangolazione_dual(const int q,const int b,const int c,PolyhedralMesh& mesh,PolyhedralMesh& meshOutput);

//carico
void Assemble_3D(PolyhedralMesh& meshTriangulated);
	
	
	
//********************* shortest

struct minimo_output {
		unsigned int num_lati;
		double length;
        minimo_output(unsigned int nEdges, double j);
	};
	
	double euclidian_dist(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2); //utilizziamo id per farlo
	MatrixXi adj_M(const PolyhedralMesh& mesh); //adjacency matrix
	
	//finder
	minimo_output minimo_F_Dijkstra(
	   PolyhedralMesh& mesh,
	   const MatrixXi& adjMatrix,
	   const unsigned int startVertexId_real,
	   const unsigned int endVertexId_real
	);
	

	
//duale**************************

void dual_calc(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);
	
Vector3d baricenter_Fetch(const PolyhedralMesh& meshTriangulated, const unsigned int faceId);
map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildMap_LF(const PolyhedralMesh& meshTriangulated);
map<unsigned int, vector<unsigned int>> buildMap_VF(const PolyhedralMesh& meshTriangolazione);
map<unsigned int, vector<unsigned int>> buildMap_VE(const PolyhedralMesh& meshTriangolazione);
	
void proj_sphere(PolyhedralMesh& meshTriangolazione);
	
	
	
	
//********************************************	
//esportazione

   void WriteCell0Ds(const PolyhedralMesh& mesh);
   void WriteCell1Ds(const PolyhedralMesh& mesh);
   void WriteCell2Ds(const PolyhedralMesh& mesh);
   void WriteCell3Ds(const PolyhedralMesh& mesh);
   
   void ExportParaview(const PolyhedralMesh& meshTriangulated);
   void MeshTriOutput(const PolyhedralMesh& meshTriangulated);


	
	
	}