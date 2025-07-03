#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;


vector<vector<int>> facce_icosaedro = {
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