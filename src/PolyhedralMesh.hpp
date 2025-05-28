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
	  // Matrice 4x3: 4 vertici, ciascuno con coordinate x, y, z
	  double Tetrahedron_Vertices[4][3] = {
        {  invSqrt3,  invSqrt3,  invSqrt3 },   // Vertice A
        {  invSqrt3, -invSqrt3, -invSqrt3 },   // Vertice B
        { -invSqrt3,  invSqrt3, -invSqrt3 },   // Vertice C
        { -invSqrt3, -invSqrt3,  invSqrt3 }    // Vertice D
		};
	
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(4, 3);
		vertices <<  invSqrt3,  invSqrt3,  invSqrt3 ,   // Vertice A
         invSqrt3, -invSqrt3, -invSqrt3 ,   // Vertice B
        -invSqrt3,  invSqrt3, -invSqrt3 ,   // Vertice C
         -invSqrt3, -invSqrt3,  invSqrt3;
		
	  int Tetrahedron_Edges[6][2] = {
		{0, 1},  // AB
		{0, 2},  // AC
		{0, 3},  // AD
		{1, 2},  // BC
		{1, 3},  // BD
		{2, 3}   // CD
		}
		
		Eigen::MatrixXi edges  = Eigen::MatrixXi::Zeros(6, 2);
		edges << 
	  int Tetrahedron_Face[4][5] = {
		{ 0, 3, 3, 0, 2, 1, {1, 3, 0} },  // Faccia A-C-B
		{ 1, 3, 3, 0, 1, 3}, {0, 4, 2} },  // Faccia A-B-D
		{ 2, 3, 3, {0, 3, 2}, {2, 5, 1} },  // Faccia A-D-C
		{ 3, 3, 3, {1, 2, 3}, {3, 5, 4} }   // Faccia B-C-D
		}
		
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
	  int Tetrahedron_Polygonal[1]= {
		  {0, 12, 6, 4, {0,1,2,3,4,5,6,7,8,9,10,11}, {0,1,2,3,4,5,6}, {0,1,2,3,4}}
	  };
	}; //end of struct
	
	struct Octahedron {
	// Matrice 6x3: 6 vertici, ciascuno con coordinate x, y, z
		double Vert_octahedron[6][3] = {
        {  1.0,  0.0,  0.0 },  // Vertice 1
        { -1.0,  0.0,  0.0 },  // Vertice 2
        {  0.0,  1.0,  0.0 },  // Vertice 3
        {  0.0, -1.0,  0.0 },  // Vertice 4
        {  0.0,  0.0,  1.0 },  // Vertice 5 
        {  0.0,  0.0, -1.0 }   // Vertice 6
		}
		
		int Octahedron_Edges[12][2] = {
		{0, 2}, {0, 3}, {0, 4}, {0, 5}, 
		{1, 2}, {1, 3}, {1, 4}, {1, 5},
		/*{2, 0}, {2, 1},*/ {2, 4}, {2, 5},
		/*{3, 0}, {3, 1},*/ {3, 4}, {3, 5},
		/*{4, 0}, {4, 1}, {4, 2}, {4, 3},*/
		/*{5, 0}, {5, 1}, {5, 2}, {5, 3},*/
		}
		
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
		
		int Octahedron_Polygonal[1]= {
		  {0, 6, 12, 8, {0,1,2,3,4,5}, {0,1,2,3,4,5,6,7,8,9,10,11}, {0,1,2,3,4,5,6,7}}
		}	
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
    double Vert_icosahedron[12][3]; //altrimenti sarebbero locali dei init
	
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
	    init_Icosahedron();
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

   void init_Icosahedron()
    {
	int index=0; // so we can initialize it in different contexts
    // Matrice 12x3 per i vertici
    //double Vert_icosahedron[12][3];

    // (0, ±1, ±φ)
    for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
            double x = 0.0;
            double y = i * 1.0;
            double z = j * phi;
            double len = sqrt(x*x + y*y + z*z);
            Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;
            index++;
        }
    }
    // (±1, ±φ, 0)
    for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
            double x = i * 1.0;
            double y = j * phi;
            double z = 0.0;
            double len = sqrt(x*x + y*y + z*z);
            Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;
            index++;
        }
    }
    // (±φ, 0, ±1)
    for (int i = -1; i <= 1; i += 2) {
        for (int j = -1; j <= 1; j += 2) {
            double x = i * phi;
            double y = 0.0;
            double z = j * 1.0;
            double len = sqrt(x*x + y*y + z*z);
            Vert_icosahedron[index][0] = x / len;
            Vert_icosahedron[index][1] = y / len;
            Vert_icosahedron[index][2] = z / len;
            index++;
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

	