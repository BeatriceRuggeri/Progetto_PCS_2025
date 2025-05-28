#pragma once

#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {
	
    struct Edge
  {
    unsigned int id;
    unsigned int origin;
    unsigned int end;
  };//end of edge struct
  
  struct Tetrahedron {
	  
	  const double invSqrt3 = 1.0 / sqrt(3.0); // 0.57735
	  /*
	  // Matrice 4x3: 4 vertici, ciascuno con coordinate x, y, z
	  double Tetrahedron_Vertices[4][3] = {
        {  invSqrt3,  invSqrt3,  invSqrt3 },   // Vertice A
        {  invSqrt3, -invSqrt3, -invSqrt3 },   // Vertice B
        { -invSqrt3,  invSqrt3, -invSqrt3 },   // Vertice C
        { -invSqrt3, -invSqrt3,  invSqrt3 }    // Vertice D
		};
	*/
	
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(4, 4);
		vertices << 0, invSqrt3, invSqrt3, invSqrt3, // Vertice A
					1, invSqrt3, -invSqrt3, -invSqrt3, // Vertice B
					2, -invSqrt3, invSqrt3, -invSqrt3, // Vertice C
					3, -invSqrt3, -invSqrt3, invSqrt3; // Vertice D 
		
	/*	
	  int Tetrahedron_Edges[6][2] = {
		{0, 1},  // AB
		{0, 2},  // AC
		{0, 3},  // AD
		{1, 2},  // BC
		{1, 3},  // BD
		{2, 3}   // CD
		}
	*/	
		Eigen::MatrixXi edges  = Eigen::MatrixXi::Zeros(6, 2);
		edges << 0, 0, 1,
				 1, 0, 2,
				 2, 0, 3,
				 3, 1, 2,
				 4, 1, 3,
				 5, 2, 3;
	/*
	  int Tetrahedron_Face[4][5] = {
		{ 0, 3, 3, 0, 2, 1, {1, 3, 0} },  // Faccia A-C-B
		{ 1, 3, 3, 0, 1, 3}, {0, 4, 2} },  // Faccia A-B-D
		{ 2, 3, 3, {0, 3, 2}, {2, 5, 1} },  // Faccia A-D-C
		{ 3, 3, 3, {1, 2, 3}, {3, 5, 4} }   // Faccia B-C-D
		}
	*/	
		//id, num_vertici, num_lati, id_vertici(3), id_lati(3)
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(4, 7);
		faces << 0, 0, 2, 1, 1, 3, 0,
				 1, 0, 1, 3, 0, 4, 2,
				 2, 0, 3, 2, 2, 5, 1,
				 3, 1, 2, 3, 3, 5, 4;
		
	/*
		int Tetrahedron_vert_face[4][3] = {
			{0, 2, 1},
			{0, 1, 3},
			{0, 3, 2},
			{1, 2, 3} }
			
		int Tetrahedron_lati_face[4][3] = {
			{1, 3, 0},
			{0, 4, 2},
			{2, 5, 1},
			{3, 5, 4}
		}
	
	*/
	//Id, num_vertici, num_lati, num_facce, id_vertici(12), id_lati(6), id_facce(4)
	  Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(28)
	  polygonal << 0, 12, 6, 4, 0,1,2,3,4,5,6,7,8,9,10,11, 0,1,2,3,4,5,6, 0,1,2,3;
	  
	  /*
	  int Tetrahedron_Polygonal[1]= {
		  {0, 12, 6, 4, {0,1,2,3,4,5,6,7,8,9,10,11}, {0,1,2,3,4,5,6}, {0,1,2,3,4}}
	  };
	  */
	}; //end of struct
	
	struct Octahedron {
		/*
	// Matrice 6x3: 6 vertici, ciascuno con coordinate x, y, z
		double Vert_octahedron[6][3] = {
        {  1.0,  0.0,  0.0 },  // Vertice 1
        { -1.0,  0.0,  0.0 },  // Vertice 2
        {  0.0,  1.0,  0.0 },  // Vertice 3
        {  0.0, -1.0,  0.0 },  // Vertice 4
        {  0.0,  0.0,  1.0 },  // Vertice 5 
        {  0.0,  0.0, -1.0 }   // Vertice 6
		}
		*/
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(6, 3);
		vertices << 0, 1.0, 0.0, 0.0,
					1, -1.0, 0.0, 0.0,
					2, 0.0, 1.0, 0.0, 
					3, 0.0, -1.0, 0.0,
					4, 0.0, 0.0, 1.0, 
					5, 0.0, 0.0, -1.0;
		/*
		int Octahedron_Edges[12][2] = {
		{0, 2}, {0, 3}, {0, 4}, {0, 5}, 
		{1, 2}, {1, 3}, {1, 4}, {1, 5},
		{2, 0}, {2, 1}, {2, 4}, {2, 5},
		{3, 0}, {3, 1}, {3, 4}, {3, 5},
		{4, 0}, {4, 1}, {4, 2}, {4, 3},
		{5, 0}, {5, 1}, {5, 2}, {5, 3},
		}
		*/
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(12, 2);
		edges << 0, 0, 2,
				 1, 0, 3, 
				 2, 0, 4, 
				 3, 0, 5, 
				 4, 1, 2, 
				 5, 1, 3, 
				 6, 1, 4, 
				 7, 1, 5,
				 8, 2, 4, 
				 9, 2, 5,
				 10, 3, 4,
				 11, 3, 5;
		/*
		int Octahedron_Face[8][5] = {
		{ 0, 3, 3, {0, 2, 4}, {0, 2, 8} },  
		{ 1, 3, 3, {2, 1, 4}, {4, 8, 6} },  
		{ 2, 3, 3, {1, 3, 4}, {5, 6, 10} },  
		{ 3, 3, 3, {3, 0, 4}, {1, 10, 2} },
		{ 4, 3, 3, {2, 0, 5}, {0, 9, 3} },
		{ 5, 3, 3, {1, 2, 5}, {4, 7, 9} },
		{ 6, 3, 3, {3, 1, 5}, {5, 11, 7} },
		{ 7, 3, 3, {0, 3, 5}, {1, 3, 11} }
		}
		*/
		
		//id, num_vertici, num_lati, id_vertici(3), id_lati(3)
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(8, 7);
		faces << 0, 0, 2, 4, 0, 2, 8,  
				 1, 2, 1, 4, 4, 8, 6,  
				 2, 1, 3, 4, 5, 6, 10,  
				 3, 3, 0, 4, 1, 10, 2,
				 4, 2, 0, 5, 0, 9, 3,
				 5, 1, 2, 5, 4, 7, 9,
				 6, 3, 1, 5, 5, 11, 7,
				 7, 0, 3, 5, 1, 3, 11;
				 
		//Id, num_vertici, num_lati, num_facce, id_vertici(6), id_lati(12), id_facce(8)
		Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(30)
		polygonal << 0, 6, 12, 8, 0,1,2,3,4,5, 0,1,2,3,4,5,6,7,8,9,10,11, 0,1,2,3,4,5,6,7;
		
		/*
		int Octahedron_Polygonal[1]= {
		  {0, 6, 12, 8, {0,1,2,3,4,5}, {0,1,2,3,4,5,6,7,8,9,10,11}, {0,1,2,3,4,5,6,7}}
		}
		*/		
	}; //end of struct
	
	struct Icosahedron {
		const double phi = (1.0 + sqrt(5.0)) / 2.0;
		const double invPhi = 1.0 / phi;
		
		void init_Icosahedron() {
			int index=0; // so we can initialize it in different contexts
			Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(12, 3); // Matrice 12x3 per i vertici
			
			// (0, ±1, ±φ)
			for (int i = -1; i <= 1; i += 2) {
				for (int j = -1; j <= 1; j += 2) {
					double x = 0.0;
					double y = i * 1.0;
					double z = j * phi;
					double len = sqrt(x*x + y*y + z*z); 
					vertices.row(index) << index, x/len, y/len, z/len;     
					index++;
			/*Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;*/
					}
			}
			// (±1, ±φ, 0)
			for (int i = -1; i <= 1; i += 2) {
				for (int j = -1; j <= 1; j += 2) {
					double x = i * 1.0;
					double y = j * phi;
					double z = 0.0;
					double len = sqrt(x*x + y*y + z*z);
					vertices.row(index)<< index, x/len, y/len, z/len;
					index++;
			
            /*Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;*/
           
				}
			}
			// (±φ, 0, ±1)
			for (int i = -1; i <= 1; i += 2) {
				for (int j = -1; j <= 1; j += 2) {
					double x = i * phi;
					double y = 0.0;
					double z = j * 1.0;
					double len = sqrt(x*x + y*y + z*z);
					vertices.row(index)<< index, x/len, y/len, z/len;
					index++;
            /*Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;*/
				}
			}
		}
		/*
		//creo i lati
		for (int i = 0; i < 12; ++i) {
			for (int j = i+1; j < 12; ++j) {
				double dist = (vertices.row(i) - vertices.row(j)).norm();
				if (dist < threshold) {
					edges.push_back({i,j});
				}
			}
		}
		
		// Stampa lati
		std::cout << "Lati (archi) dell'icosaedro (totali: " << edges.size() << "):\n";
		for (auto &e : edges) {
			std::cout << e.first << " - " << e.second << "\n";
		}*/

		// Definizione dei 30 lati (archi)
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(30, 2);
		edges << 0, 0, 1, //AB 
				 1, 0, 4, //AE 
				 2, 0, 6, //AG 
				 3, 0, 8, //AI 
				 4, 0, 10, //AM 
				 5, 1, 5, //BF 
				 6, 1, 8, //BI 
				 7, 1, 7, //BH 
				 8, 1, 10, //BM 
				 9, 2, 3, //CD 
				 10, 2, 4, //CE 
				 11, 2, 6, //CG 
				 12, 2, 9, //CL 
				 13, 2, 11, //CN 
				 14, 3, 5, //DF 
				 15, 3, 7, //DH 
				 16, 3, 9, //DL 
				 17, 3, 11, //DN 
				 18, 4, 5, //EF 
				 19, 4, 6, //EG 
				 20, 4, 8, //EI 
				 21, 4, 9, //EL 
				 22, 5, 7, //FH 
				 23, 5, 8, //FI 
				 24, 5, 9, //FL 
				 25, 6, 7, //GH 
				 26, 6, 10, //GM 
				 27, 6, 11, //GN 
				 28, 8, 9, //IL 
				 29, 10, 11; //MN 
				 
		//Definisco le faccie  (togli num_vertici e num_lati)
		//id, num_vertici, num_lati, id_vertici(3), id_lati(3)
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(20, 9);
		faces << 0, 3, 3, ,
				 1, 3, 3, 
				 2, 3, 3,
				 3, 3, 3,
				 4, 3, 3,
				 5, 3, 3,
				 6, 3, 3,
				 7, 3, 3,
				 8, 3, 3,
				 9, 3, 3,
				 10, 3, 3,
				 11, 3, 3,
				 12, 3, 3,
				 13, 3, 3,
				 14, 3, 3,
				 15, 3, 3,
				 16, 3, 3,
				 17, 3, 3,
				 18, 3, 3,
				 19, 3, 3,
		

    std::vector<std::array<int,3>> face_vertices = {
        {}, {0, 5, 10}, {2, 4, 9}, {2, 11, 5}, {1, 6, 8},
        {1, 10, 7}, {3, 9, 6}, {3, 7, 11}, {0, 10, 8}, {1, 8, 10},
        {2, 9, 11}, {3, 11, 9}, {4, 8, 6}, {5, 7, 10}, {6, 7, 11},
        {7, 9, 11}, {4, 5, 9}, {5, 6, 7}, {6, 4, 7}, {8, 10, 6}
    };
		
	}; //end of struct


    struct PolyhedralMesh
{
    const double invSqrt3 = 1.0 / sqrt(3.0); // 0.57735

    // Matrice 8x3: 8 vertici, ciascuno con coordinate x, y, z
    double Vert_cube[8][3] = {
        {  invSqrt3,  invSqrt3,  invSqrt3 },  // Vertice 1
        {  invSqrt3,  invSqrt3, -invSqrt3 },  // Vertice 2
        {  invSqrt3, -invSqrt3,  invSqrt3 },  // Vertice 3
        {  invSqrt3, -invSqrt3, -invSqrt3 },  // Vertice 4
        { -invSqrt3,  invSqrt3,  invSqrt3 },  // Vertice 5
        { -invSqrt3,  invSqrt3, -invSqrt3 },  // Vertice 6
        { -invSqrt3, -invSqrt3,  invSqrt3 },  // Vertice 7
        { -invSqrt3, -invSqrt3, -invSqrt3 }   // Vertice 8
    }

    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    const double invPhi = 1.0 / phi;
    
    double Vert_dodecahedron[20][3];
    //double Vert_icosahedron[12][3]; //altrimenti sarebbero locali dei init
	
	int cube_edges[12][2] = {
		{0,1}, {0,2}, {0,4},
        {1,3}, {1,5},
        {2,3}, {2,6},
        {3,7},
        {4,5}, {4,6},
        {5,7},
        {6,7}
	   }
	   
    //constructor
    PolyhedralMesh(){
	    //init_Icosahedron();
	    init_Dodecahedron();

	    
    }
	
	
	
    

// in teoria possiamo mettere questo su un PolyhedralMesh.cpp eventualmente boh

//************************* initialization ********************//


    // Matrice 20x3: ogni riga è un vertice (x, y, z)
    //double Vert_dodecahedron[20][3];
    //int index = 0; //risking out of bounds array writing if we use it for different functions

    // Vertici (±1, ±1, ±1) — 8 vertici
   void init_Dodecahedron()
    {
	int index=0;
    for (int sx = -1; sx <= 1; sx += 2)
        for (int sy = -1; sy <= 1; sy += 2)
            for (int sz = -1; sz <= 1; sz += 2) {
                double x = sx * 1.0;
                double y = sy * 1.0;
                double z = sz * 1.0;
                double len = std::sqrt(x*x + y*y + z*z); //normalization
                Vert_dodecahedron[index][0] = x / len;
                Vert_dodecahedron[index][1] = y / len;
                Vert_dodecahedron[index][2] = z / len;
                index++;
            }

    // Vertici (0, ±1/φ, ±φ) + permutazioni — 12 vertici
    double coords[3][2] = {
        {0.0, 0.0},
        {invPhi, -invPhi},
        {phi, -phi}
    };

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                for (int l = 0; l < 2; ++l) {
                    double x = coords[i][j];
                    double y = coords[(i + 1) % 3][k];
                    double z = coords[(i + 2) % 3][l];
                    double len = std::sqrt(x*x + y*y + z*z); // normalization 
                    Vert_dodecahedron[index][0] = x / len;
                    Vert_dodecahedron[index][1] = y / len;
                    Vert_dodecahedron[index][2] = z / len;
                    index++;
                 }
             }
         }
      }
    }

   
   
  /* void Cube_Edges(){
	   
	  
	   
	   vector<pair<int, int>> cube_edges = {
		{0,1}, {0,2}, {0,4},
        {1,3}, {1,5},
        {2,3}, {2,6},
        {3,7},
        {4,5}, {4,6},
        {5,7},
        {6,7}
	   };
	   for (std::size_t i=0; i< cube_edges.size(); ++i) {
		   Cell1Ds.push_back({i, cube_edges[i].first, cube_edges[i].second});
	   }

    }
	
	*/
	
	
 }; // end of PolyhedralMesh struct
 
}//end library

	