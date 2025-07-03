#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;

vector<vector<int>> spigoli_icosaedro = {
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