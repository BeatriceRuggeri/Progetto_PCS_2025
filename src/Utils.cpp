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
        if (!ExportCell0Ds(mesh))
            return false;

        if (!ExportCell1Ds(mesh))
            return false;

        if (!ExportCell2Ds(mesh))
            return false;
            
        if (!ExportCell3DS(mesh))
            return false;

        return true;
    }
    
    
    
    
    
    
    
//*********************************************** Importazione ********************* : Ancora d'adattare 

//**********CELL 0DS************

  
    bool ImportCell0Ds(PolyhedrallMesh& mesh)
    {
        ifstream file("./Cell0Ds.txt");

        if (file.fail())
            return false;

        list<string> listLines;
        string line;
        while (getline(file, line))
            listLines.push_back(line);

        file.close();

        // remove header
        listLines.pop_front();

        mesh.NumCell0Ds = listLines.size();

        if (mesh.NumCell0Ds == 0)
        {
            cerr << "There is no cell 0D" << endl;
            return false;
        }

        mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
        mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);  // Correct size ?
 
        for (const string& line : listLines)
        {
          istringstream converter(line);
          unsigned int id;
          unsigned int marker;
          Vector2d coord;

          // Read the values
          char tmp;
          converter >>  id >>tmp>> marker >>tmp>>mesh.Cell0DsCoordinates(0,id)>>tmp>>mesh.Cell0DsCoordinates(1,id); 
         mesh.Cell0DsId.push_back(id);

         // markers
         if (marker != 0)
          {
            const auto it = mesh.MarkerCell0Ds.find(marker);
            if (it == mesh.MarkerCell0Ds.end())
              {
                mesh.MarkerCell0Ds.insert({marker, {id}});
              }
         else
            {
            it->second.push_back(id);
          }
       }
}


        return true;
    }

// *********************************** CELL 1 DS ***************************

    bool ImportCell1Ds(PolygonalMesh& mesh)
    {
        ifstream file("./Cell1Ds.txt");

        if (file.fail())
            return false;

        list<string> listLines;
        string line;
        while (getline(file, line))
            listLines.push_back(line);

        file.close();

        // remove header
        listLines.pop_front();

        mesh.NumCell1Ds = listLines.size();

        if (mesh.NumCell1Ds == 0)
        {
            cerr << "There is no cell 1D" << endl;
            return false;
        }

        mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
        mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

       
        for (const string& line : listLines)
        {
            istringstream converter(line);
            unsigned int id;
            unsigned int marker;
            Vector2i vertices;
            char tmp;

            converter >>  id >>tmp>> marker >>tmp >> mesh.Cell1DsExtrema(0, id) >> tmp >>mesh.Cell1DsExtrema(1, id);
            mesh.Cell1DsId.push_back(id);
            
            if (marker != 0)
            {
                const auto it = mesh.MarkerCell1Ds.find(marker);
                if (it == mesh.MarkerCell1Ds.end())
                {
                    mesh.MarkerCell1Ds.insert({marker, {id}});
                }
                else
                {
                    it->second.push_back(id);
                }
            }
           
        }

        return true;
    }

//************************************************* CELL 2DS********************************************

bool ImportCell2Ds(PolyhedrallMesh& mesh){
	
	
	
	
	
	return true;
}

//************************************************* CELL 3DS********************************************

bool ImportCell3Ds(PolyhedrallMesh& mesh){
	
	
	
	return true;
}