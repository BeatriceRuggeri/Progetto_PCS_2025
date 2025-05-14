#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
using namespace PolyhedrallLibrary;





namespace PolyhedralLibrary
{
    // ***************************************************************************
    bool ImportMesh(PolygonalMesh& mesh)
    {
        if (!ImportCell0Ds(mesh))
            return false;

        if (!ImportCell1Ds(mesh))
            return false;

        if (!ImportCell2Ds(mesh))
            return false;

        return true;
    }
