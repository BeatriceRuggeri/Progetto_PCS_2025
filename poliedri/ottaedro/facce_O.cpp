#include <iostream>
#include <cmath>
#include "Eigen/Eigen"

using namespace std;


vector<vector<int>> facce_ottaedro = {
    {0, 0, 2, 4, 0, 8, 2}, //A-C-E -- AC,CE,AE
    {1, 0, 3, 4, 2, 1, 10}, //A-D-E -- AE,AD,DE 
    {2, 1, 3, 4, 10, 5, 6}, //B-D-E -- DE,BD,BE
    {3, 1, 2, 4, 6, 8, 4}, //B-C-E -- BE,CE,BC
    {4, 1, 2, 5, 4, 7, 9}, //B-C-F -- BC,BF,CF
    {5, 2, 0, 5, 9, 0, 3}, //A-C-F -- CF,AC,AF
    {6, 0, 3, 5, 3, 1, 11}, //A-D-F -- AF,AD,DF
    {7, 3, 1, 5, 11, 5, 7} //B-D-F -- DF,BD,BF
};