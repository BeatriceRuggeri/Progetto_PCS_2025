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
	
		vector<vector<double>> vertices = {
			{0, invSqrt3, invSqrt3, invSqrt3}, // A
			{1, invSqrt3, -invSqrt3, -invSqrt3}, // B
			{2, -invSqrt3, invSqrt3, -invSqrt3}, // C
			{3, -invSqrt3, -invSqrt3, invSqrt3} // D
		};
		/*
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(4, 4);
		vertices << 0, invSqrt3, invSqrt3, invSqrt3, // A
					1, invSqrt3, -invSqrt3, -invSqrt3, // B
					2, -invSqrt3, invSqrt3, -invSqrt3, // C
					3, -invSqrt3, -invSqrt3, invSqrt3; // D 
		*/
	
		vector<vector<int>> edges = {
			{0, 0, 1}, // AB
			{1, 0, 2}, // AC
			{2, 0, 3}, // AD
			{3, 1, 2}, // BC
			{4, 1, 3}, // BD
			{5, 2, 3} // CD
		};
		/*
		Eigen::MatrixXi edges  = Eigen::MatrixXi::Zeros(6, 3);
		edges << 0, 0, 1, // AB
				 1, 0, 2, // AC
				 2, 0, 3, // AD
				 3, 1, 2, // BC
				 4, 1, 3, // BD
				 5, 2, 3; // CD
		*/
	
		vector<vector<int>> faces = {
			{0, 0, 2, 1, 1, 3, 0}, //A-B-C -- AC,BC,AB
			{1, 0, 1, 3, 0, 4, 2}, //A-B-D -- AB,BD,AD
			{2, 0, 3, 2, 2, 1, 5}, //A-C-D -- AD,AC,CD
			{3, 1, 2, 3, 5, 3, 4} //B-C-D -- CD,BC,BD
		};
			
		/*
		//id, num_vertici, num_lati, id_vertici(3), id_lati(3)
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(4, 7);
		faces << 0, 0, 2, 1, 1, 3, 0, //A-B-C -- AC,BC,AB
				 1, 0, 1, 3, 0, 4, 2, //A-B-D -- AB,BD,AD
				 2, 0, 3, 2, 2, 1, 5, //A-C-D -- AD,AC,CD
				 3, 1, 2, 3, 5, 3, 4; //B-C-D -- CD,BC,BD */
				 
		
	//Id, num_vertici, num_lati, num_facce, id_vertici(12), id_lati(6), id_facce(4)
	  Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(28)
	  polygonal << 0, 12, 6, 4, 0,1,2,3,4,5,6,7,8,9,10,11, 0,1,2,3,4,5,6, 0,1,2,3;

	}; //end of struct
	
	struct Octahedron {
	
		vector<vector<double>> vertices = {
			{0, 1.0, 0.0, 0.0}, // A
			{1, -1.0, 0.0, 0.0}, // B
			{2, 0.0, 1.0, 0.0}, // C
			{3, 0.0, -1.0, 0.0}, // D
			{4, 0.0, 0.0, 1.0}, // E
			{5, 0.0, 0.0, -1.0} // F
		};
		/*
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(6, 4);
		vertices << 0, 1.0, 0.0, 0.0, // A
					1, -1.0, 0.0, 0.0, // B
					2, 0.0, 1.0, 0.0, // C
					3, 0.0, -1.0, 0.0, // D
					4, 0.0, 0.0, 1.0, // E
					5, 0.0, 0.0, -1.0; // F */

		vector<vector<int>> edges = {
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
		
		/*	
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(12, 3);
		edges << 0, 0, 2, // AC
				 1, 0, 3, // AD
				 2, 0, 4, // AE
				 3, 0, 5, // AF
				 4, 1, 2, // BC
				 5, 1, 3, // BD
				 6, 1, 4, // BE
				 7, 1, 5, // BF
				 8, 2, 4, // CE
				 9, 2, 5, // CF
				 10, 3, 4, // DE
				 11, 3, 5; // DF */
		
		vector<vector<int>> faces = {
			{0, 0, 2, 4, 0, 8, 2}, //A-C-E -- AC,CE,AE
			{1, 0, 3, 4, 2, 1, 10}, //A-D-E -- AE,AD,DE 
			{2, 1, 3, 4, 10, 5, 6}, //B-D-E -- DE,BD,BE
			{3, 1, 2, 4, 6, 8, 4}, //B-C-E -- BE,CE,BC
			{4, 1, 2, 5, 4, 7, 9}, //B-C-F -- BC,BF,CF
			{5, 2, 0, 5, 9, 0, 3}, //A-C-F -- CF,AC,AF
			{6, 0, 3, 5, 3, 1, 11}, //A-D-F -- AF,AD,DF
			{7, 3, 1, 5, 11, 5, 7} //B-D-F -- DF,BD,BF
		};
		/*
		//id, id_vertici(3), id_lati(3)
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(8, 7);
		faces << 0, 0, 2, 4, 0, 8, 2, //A-C-E -- AC,CE,AE
				 1, 0, 3, 4, 2, 1, 10, //A-D-E -- AE,AD,DE 
				 2, 1, 3, 4, 10, 5, 6, //B-D-E -- DE,BD,BE
				 3, 1, 2, 4, 6, 8, 4, //B-C-E -- BE,CE,BC
				 4, 1, 2, 5, 4, 7, 9, //B-C-F -- BC,BF,CF
				 5, 2, 0, 5, 9, 0, 3, //A-C-F -- CF,AC,AF
				 6, 0, 3, 5, 3, 1, 11, //A-D-F -- AF,AD,DF
				 7, 3, 1, 5, 11, 5, 7; //B-D-F -- DF,BD,BF */
		
		//Id, num_vertici, num_lati, num_facce, id_vertici(6), id_lati(12), id_facce(8)
		Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(30)
		polygonal << 0, 6, 12, 8, 0,1,2,3,4,5, 0,1,2,3,4,5,6,7,8,9,10,11, 0,1,2,3,4,5,6,7;
			
	}; //end of struct
	
	struct Icosahedron {
		const double Phi = (1.0 + sqrt(5.0)) / 2.0;
		const double Raz_Phi = sqrt(Phi^2 + 1);
		const double OP = 1.0/Raz_Phi; //OP=one su Phi
		const double GP = Phi/Raz_Phi; // GP = giusto Phi= Phi/Raz_Phi
		
		vector<vector<double>> vertices = {
			{0, 0.0, OP, GP}, // A
			{1, 0.0, -OP, GP}, // B
			{2, 0.0, OP, -GP}, // C
			{3, 0.0, -OP, -GP}, // D
			{4, OP, GP, 0.0}, // E
			{5, OP, -GP, 0.0}, // F
			{6, -OP, GP, 0.0}, // G
			{7, -OP, -GP, 0.0}, // H
			{8, GP, 0.0, OP}, // I
			{9, GP, 0.0, -OP}, //L
			{10, -GP, 0.0, OP}, // M
			{11, -GP, 0.0, -OP} //N
		};
		/*
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(12, 4);
		vertices << 0, 0.0, OP, GP, // A
					1, 0.0, -OP, GP, // B
					2, 0.0, OP, -GP, // C
					3, 0.0, -OP, -GP, // D
					4, OP, GP, 0.0, // E
					5, OP, -GP, 0.0, // F
					6, -OP, GP, 0.0, // G
					7, -OP, -GP, 0.0, // H
					8, GP, 0.0, OP, // I
					9, GP, 0.0, -OP, //L
					10, -GP, 0.0, OP, // M
					11, -GP, 0.0, -OP; //N */

		vector<vector<int>> edges = {
			{0, 0, 1}, //AB 
			{1, 0, 4}, //AE 
			{2, 0, 6}, //AG 
			{3, 0, 8}, //AI 
			{4, 0, 10}, //AM 
			{5, 1, 5}, //BF 
			{6, 1, 8}, //BI 
			{7, 1, 7}, //BH 
			{8, 1, 10}, //BM 
			{9, 2, 3}, //CD 
			{10, 2, 4}, //CE 
			{11, 2, 6}, //CG 
			{12, 2, 9}, //CL 
			{13, 2, 11}, //CN 
			{14, 3, 5}, //DF 
			{15, 3, 7}, //DH 
			{16, 3, 9}, //DL 
			{17, 3, 11}, //DN  
			{18, 4, 6}, //EG 
			{19, 4, 8}, //EI 
			{20, 4, 9}, //EL 
			{21, 5, 7}, //FH 
			{22, 5, 8}, //FI 
			{23, 5, 9}, //FL 
			{24, 6, 10}, //GM 
			{25, 6, 11}, //GN 
			{26, 7, 10}, //HM
			{27, 7, 11}, //HN
			{28, 8, 9}, //IL 
			{29, 10, 11} //MN
		};
		
		/*
		// Definizione dei 30 lati (archi)
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(30, 3);
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
				 18, 4, 6, //EG 
				 19, 4, 8, //EI 
				 20, 4, 9, //EL 
				 21, 5, 7, //FH 
				 22, 5, 8, //FI 
				 23, 5, 9, //FL 
				 24, 6, 10, //GM 
				 25, 6, 11, //GN 
				 26, 7, 10, //HM
				 27, 7, 11, //HN
				 28, 8, 9, //IL 
				 29, 10, 11; //MN */		 
				 
		//Definisco le faccie  
		//id, id_vertici(3), id_lati(3)
		vector<vector<int>> faces = {
			{0, 0, 1, 8, 3, 6, 0}, // A-B-I - AI, BI, AB                
			{1, 0, 1, 10, 0, 8, 4}, // A-B-M - AB, BM, AM,
			{2, 0, 6, 10, 4, 26, 2}, // A-G-M - AM, GM, AG
			{3, 0, 4, 6, 2, 19, 1}, // A-E-G - AG, EG, AE, 
			{4, 0, 4, 8, 1, 3, 19}, // A-I-E - AE, AI, EI 	
			{5, 4, 8, 9, 19, 20, 28}, // E-I-L - EI, EL, IL,
			{6, 5, 8, 9, 28, 23, 22}, // F-I-L - IL, FL, FI, 
			{7, 1, 5, 8, 23, 6, 5}, // B-I-F - IF, BI, BF, 
			{8, 1, 5, 7, 5, 21, 7}, // B-H-F - BF, FH, BH, 
			{9, 1, 7, 10, 7, 8, 26}, // B-H-M - BH, BM, HM, 
			{10, 7, 10, 11, 26, 27, 29}, // H-M-N - HM, HN, MN, 
			{11, 6, 10, 11, 29, 24, 25}, // G-M-N - MN, GM, GN
			{12, 2, 6, 11, 25, 13, 11},  // C-N-G - GN, CN, CG,
			{13, 2, 4, 6, 11, 18, 10}, // C-E-G - CG, EG, CE,	
			{14, 2, 4, 9, 10, 20, 12}, // C-E-L - CE, EL, CL,
			{15, 2, 3, 9, 12, 9, 16}, // C-D-L - CL, CD, DL,
			{16, 3, 5, 9, 16, 23, 14}, // F-L-D - DL, FL, DF,
			{17, 3, 5, 7, 14, 21, 15}, // D-F-H - DF, FH, DH,
			{18, 3, 7, 11, 15, 27, 17}, // D-H-N - DH, HN, DN,
			{19, 2, 3, 11, 17, 13, 9} // C-D-N - DN, CN, CD
		};
		
		/*
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(20, 7);
		faces << 0, 0, 1, 8, 3, 6, 0, // A-B-I - AI, BI, AB                
				 1, 0, 1, 10, 0, 8, 4, // A-B-M - AB, BM, AM,
				 2, 0, 6, 10, 4, 26, 2, // A-G-M - AM, GM, AG
				 3, 0, 4, 6, 2, 19, 1, // A-E-G - AG, EG, AE, 
				 4, 0, 4, 8, 1, 3, 19, // A-I-E - AE, AI, EI 	
				 5, 4, 8, 9, 19, 20, 28, // E-I-L - EI, EL, IL,
				 6, 5, 8, 9, 28, 23, 22, // F-I-L - IL, FL, FI, 
				 7, 1, 5, 8, 23, 6, 5, // B-I-F - IF, BI, BF, 
				 8, 1, 5, 7, 5, 21, 7, // B-H-F - BF, FH, BH, 
				 9, 1, 7, 10, 7, 8, 26, // B-H-M - BH, BM, HM, 
				 10, 7, 10, 11, 26, 27, 29, // H-M-N - HM, HN, MN, 
				 11, 6, 10, 11, 29, 24, 25, // G-M-N - MN, GM, GN
				 12, 2, 6, 11, 25, 13, 11,  // C-N-G - GN, CN, CG,
				 13, 2, 4, 6, 11, 18, 10, // C-E-G - CG, EG, CE,	
				 14, 2, 4, 9, 10, 20, 12, // C-E-L - CE, EL, CL,
				 15, 2, 3, 9, 12, 9, 16, // C-D-L - CL, CD, DL,
				 16, 3, 5, 9, 16, 23, 14, // F-L-D - DL, FL, DF,
				 17, 3, 5, 7, 14, 21, 15, // D-F-H - DF, FH, DH,
				 18, 3, 7, 11, 15, 27, 17, // D-H-N - DH, HN, DN,
				 19, 2, 3, 11, 17, 13, 9; // C-D-N - DN, CN, CD, */
		
		//Id, num_vertici, num_lati, num_facce, id_vertici(12), id_lati(30), id_facce(20)
		Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(66)
		polygonal << 0, 12, 30, 20, 0,1,2,3,4,5,6,7,8,9,10,11, <<
				  << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29, << 
				  << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19;
	}; //end of struct
	
	struct Cube {
		const double invSqrt3 = 1.0 / sqrt(3.0); // 0.57735
		
		vector<vector<double>> vertices = {
			{0, invSqrt3, invSqrt3, invSqrt3}, // A
			{1, invSqrt3,  invSqrt3, -invSqrt3}, // B
			{2, invSqrt3, -invSqrt3, invSqrt3}, // C
			{3, invSqrt3, -invSqrt3, -invSqrt3}, // D
			{4, -invSqrt3, invSqrt3, invSqrt3}, // E
			{5, -invSqrt3, invSqrt3, -invSqrt3}, // F
			{6, -invSqrt3, -invSqrt3, invSqrt3}, // G
			{7, -invSqrt3, -invSqrt3, -invSqrt3} // H
		};
		
		/*
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(8, 4);
		vertices << 0, invSqrt3, invSqrt3, invSqrt3, // A
					1, invSqrt3,  invSqrt3, -invSqrt3, // B
					2, invSqrt3, -invSqrt3, invSqrt3, // C
					3, invSqrt3, -invSqrt3, -invSqrt3, // D
					4, -invSqrt3, invSqrt3, invSqrt3, // E
					5, -invSqrt3, invSqrt3, -invSqrt3, // F
					6, -invSqrt3, -invSqrt3, invSqrt3, // G
					7, -invSqrt3, -invSqrt3, -invSqrt3; // H */
		
		vector<vector<int>> edges = {
			{0, 0, 1}, // AB 
			{1, 0, 2}, // AC
			{2, 0, 4}, // AE
			{3, 1, 3}, // BD  
			{4, 1, 5}, // BF  
			{5, 2, 3}, // CD
			{6, 2, 6}, // CG
			{7, 3, 7}, // DH 
			{8, 4, 5}, // EF
			{9, 4, 6}, // EG 
			{10, 5, 7}, // FH 
			{11, 6, 7} // GH
		};
		
		/*
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(12, 3);
		edges << 0, 0, 1, // AB 
				 1, 0, 2, // AC
				 2, 0, 4, // AE
				 3, 1, 3, // BD  
				 4, 1, 5, // BF  
				 5, 2, 3, // CD
				 6, 2, 6, // CG
				 7, 3, 7, // DH 
				 8, 4, 5, // EF
				 9, 4, 6, // EG 
				 10, 5, 7, // FH 
				 11, 6, 7; // GH */
		
		//id, id_vertici(4), id_lati(4)
		vector<vector<int>> faces = {
			{0, 0, 2, 4, 6, 1, 2, 6, 9}, //A-C-E-G -- AC,AE,CG,EG
			{1, 4, 5, 6, 7, 9, 11, 10, 8}, //E-F-G-H -- EG,GH,FH,EF
			{2, 0, 1, 4, 5, 8, 2, 4, 0}, //A-B-E-F -- EF,AE,BF,AB
			{3, 0, 1, 2, 3, 0, 1, 3, 5}, //A-B-C-D -- AB,AC,BD,CD
			{4, 2, 3, 6, 7, 5, 6, 11, 7}, //C-D-G-H -- CD,CG,GH,DH
			{5, 1, 3, 5, 7, 7, 3, 4, 10} //B-D-F-H -- DH,BD,BF,FH
		};
		/*
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(6, 9);
		faces << 0, 0, 2, 4, 6, 1, 2, 6, 9, //A-C-E-G -- AC,AE,CG,EG
				 1, 4, 5, 6, 7, 9, 11, 10, 8, //E-F-G-H -- EG,GH,FH,EF
				 2, 0, 1, 4, 5, 8, 2, 4, 0, //A-B-E-F -- EF,AE,BF,AB
				 3, 0, 1, 2, 3, 0, 1, 3, 5, //A-B-C-D -- AB,AC,BD,CD
				 4, 2, 3, 6, 7, 5, 6, 11, 7, //C-D-G-H -- CD,CG,GH,DH
				 5, 1, 3, 5, 7, 7, 3, 4, 10; //B-D-F-H -- DH,BD,BF,FH */
		
		//Id, num_vertici, num_lati, num_facce, id_vertici(8), id_lati(12), id_facce(6)
		Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(30)
		polygonal << 0, 8, 12, 6, 0,1,2,3,4,5,6,7, <<
				  << 0,1,2,3,4,5,6,7,8,9,10,11 << 0,1,2,3,4,5; 
		
		
	}; //end of struct
	
	struct Dodecahedron {
		
		vector<vector<double>> vertices = {
			{0, 0.6, 0, 0.8}, //A
			{1, 0.18, 0.57, 0.8}, //B
			{2, -0.48, 0.35, 0.8}, //C
			{3, -0.48, -0.35, 0.8}, //D
			{4, 0.18, -0.35, 0.8}, //E
			{5, 0.98, 0, 0.18}, //F
			{6, 0.30, 0.93, 0.18}, //G
			{7, -0.78, 0.57, 0.18}, //H
			{8, -0.78, -0.57, 0.18}, //I
			{9, 0.30, -0.93, 0.18}, //L
			{10, -0.98, 0, -0.18}, //M
			{11, -0.30, -0.93, -0.18}, //N
			{12, 0.78, -0.57, -0.18}, //O
			{13, 0.78, 0.57, -0.18}, //P
			{14, -0.30, 0.93, -0.18}, //Q
			{15, -0.6, 0, -0.8}, //R
			{16, -0.18, -0.57, -0.8}, //S
			{17, 0.48, -0.35, -0.8}, //T
			{18, 0.48, 0.35, -0.8}, //U
			{19, -0.18, 0.57, -0.8} //V
		};
		/*
		Eigen::MatrixXd vertices = Eigen::MatrixXd::Zeros(20, 4);
		vertices << 0, 0.6, 0, 0.8, //A
					1, 0.18, 0.57, 0.8, //B
					2, -0.48, 0.35, 0.8, //C
					3, -0.48, -0.35, 0.8, //D
					4, 0.18, -0.35, 0.8, //E
					5, 0.98, 0, 0.18, //F
					6, 0.30, 0.93, 0.18, //G
					7, -0.78, 0.57, 0.18, //H
					8, -0.78, -0.57, 0.18, //I
					9, 0.30, -0.93, 0.18, //L
					10, -0.98, 0, -0.18, //M
					11, -0.30, -0.93, -0.18, //N
					12, 0.78, -0.57, -0.18, //O
					13, 0.78, 0.57, -0.18, //P
					14, -0.30, 0.93, -0.18, //Q
					15, -0.6, 0, -0.8, //R
					16, -0.18, -0.57, -0.8, //S
					17, 0.48, -0.35, -0.8, //T
					18, 0.48, 0.35, -0.8, //U
					19, -0.18, 0.57, -0.8; //V	*/
		
		vector<vector<int>> edges = {
			{0, 0, 1}, //AB
			{1, 0, 4}, //AE
			{2, 0, 5}, //AF
			{3, 1, 2}, //BC
			{4, 1, 6}, //BG
			{5, 2, 3}, //CD
			{6, 2, 7}, //CH
			{7, 3, 4}, //DE
			{8, 3, 8}, //DI
			{9, 4, 9}, //EL
			{10, 5, 12}, //FO
			{11, 5, 13}, //FP
			{12, 6, 13}, //GP
			{13, 6, 14}, //GQ
			{14, 7, 10}, //HM
			{15, 7, 14}, //HQ
			{16, 8, 10}, //IM
			{17, 8, 11}, //IN
			{18, 9, 11}, //LN
			{19, 9, 12}, //LO
			{20, 10, 15}, //MR
			{21, 11, 16}, //NS
			{22, 12, 17}, //OT
			{23, 13, 18}, //PU
			{24, 14, 19}, //QV
			{25, 15, 16}, //RS
			{26, 15, 19}, //RV
			{27, 16, 17}, //ST
			{28, 17, 18}, //TU
			{29, 18, 19} //UV
		};
		
		/*
		Eigen::MatrixXi edges = Eigen::MatrixXi::Zeros(30, 3);
		edges << 0, 0, 1, //AB
				 1, 0, 4, //AE
				 2, 0, 5, //AF
				 3, 1, 2, //BC
				 4, 1, 6, //BG
				 5, 2, 3, //CD
				 6, 2, 7, //CH
				 7, 3, 4, //DE
				 8, 3, 8, //DI
				 9, 4, 9, //EL
				 10, 5, 12, //FO
				 11, 5, 13, //FP
				 12, 6, 13, //GP
				 13, 6, 14, //GQ
				 14, 7, 10, //HM
				 15, 7, 14, //HQ
				 16, 8, 10, //IM
				 17, 8, 11, //IN
				 18, 9, 11, //LN
				 19, 9, 12, //LO
				 20, 10, 15, //MR
				 21, 11, 16, //NS
				 22, 12, 17, //OT
				 23, 13, 18, //PU
				 24, 14, 19, //QV
				 25, 15, 16, //RS
				 26, 15, 19, //RV
				 27, 16, 17, //ST
				 28, 17, 18, //TU
				 29, 18, 19; //UV	*/
		
		vector<vector<int>> faces = {
			{0, 0, 1, 2, 3, 4, 0, 1, 3, 7, 5}, //A-B-C-D-E -- AB,AE,BC,DE,CD
			{1, 2, 3, 7, 8, 10, 5, 6, 14, 16, 8}, //C-D-H-I-M --CD,CH,HM,IM,DI
			{2, 3, 4, 8, 9, 11, 8, 7, 17, 18, 9}, //D-E-I-L-N -- DI,DE,IN,LN,EL
			{3, 0, 4, 5, 9, 12, 9, 1, 19, 10, 2}, //A-E-F-L-O -- EL,AE,LO,FO,AF
			{4, 0, 1, 5, 6, 13, 2, 0, 11, 12, 4}, //A-B-F-G-P -- AF,AB,FP,GP,BG
			{5, 1, 2, 6, 7, 14, 4, 3, 13, 6, 15}, //B-C-G-H-Q -- BG,BC,GQ,CH,HQ
			{6, 7, 10, 14, 15, 19, 15, 14, 24, 26, 20}, //H-M-Q-R-V -- HQ,HM,QV,RV,MR
			{7, 8, 10, 11, 15, 16, 20, 16, 17, 25, 21}, //I-M-N-R-S -- MR,IM,IN,RS,NS,
			{8, 9, 11, 12, 16, 17, 21, 18, 19, 27, 22}, //L-N-O-S-T -- NS,LN,LO,ST,OT,
			{9, 5, 12, 13, 17, 18, 22, 10, 11, 28, 23}, //F-O-P-T-U -- OT,FO,FP,TU,PU,
			{10, 6, 13, 14, 18, 19, 23, 12, 13, 24, 29}, //G-P-Q-U-V -- PU,GP,GQ,QV,UV
			{11, 15, 16, 17, 18, 19, 29, 28, 26, 25, 27} //R-S-T-U-V -- UV,UT,RV,RS,ST
		};
		/*
		Eigen::MatrixXi faces  = Eigen::MatrixXi::Zeros(12, 11);
		faces << 0, 0, 1, 2, 3, 4, 0, 1, 3, 7, 5, //A-B-C-D-E -- AB,AE,BC,DE,CD
				 1, 2, 3, 7, 8, 10, 5, 6, 14, 16, 8, //C-D-H-I-M --CD,CH,HM,IM,DI
				 2, 3, 4, 8, 9, 11, 8, 7, 17, 18, 9, //D-E-I-L-N -- DI,DE,IN,LN,EL
				 3, 0, 4, 5, 9, 12, 9, 1, 19, 10, 2, //A-E-F-L-O -- EL,AE,LO,FO,AF
				 4, 0, 1, 5, 6, 13, 2, 0, 11, 12, 4, //A-B-F-G-P -- AF,AB,FP,GP,BG
				 5, 1, 2, 6, 7, 14, 4, 3, 13, 6, 15, //B-C-G-H-Q -- BG,BC,GQ,CH,HQ
				 6, 7, 10, 14, 15, 19, 15, 14, 24, 26, 20, //H-M-Q-R-V -- HQ,HM,QV,RV,MR
				 7, 8, 10, 11, 15, 16, 20, 16, 17, 25, 21, //I-M-N-R-S -- MR,IM,IN,RS,NS,
				 8, 9, 11, 12, 16, 17, 21, 18, 19, 27, 22, //L-N-O-S-T -- NS,LN,LO,ST,OT,
				 9, 5, 12, 13, 17, 18, 22, 10, 11, 28, 23, //F-O-P-T-U -- OT,FO,FP,TU,PU,
				 10, 6, 13, 14, 18, 19, 23, 12, 13, 24, 29, //G-P-Q-U-V -- PU,GP,GQ,QV,UV
				 11, 15, 16, 17, 18, 19, 29, 28, 26, 25, 27; //R-S-T-U-V -- UV,UT,RV,RS,ST */
				 
		Eigen::VectorXi polygonal = Eigen::VectorXi::Zeros(66)
		polygonal << 0, 20, 30, 12, <<
				  << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19, <<
				  << 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29, <<
				  << 0,1,2,3,4,5,6,7,8,9,10,11;
		
		
	}; //end of struct



    struct PolyhedralMesh
{
   
::
	
 }; // end of PolyhedralMesh struct
 
}//end library

	