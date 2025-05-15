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
    
    
    
    
    
    
    
//*********************************************** Esportazione ********************* : Ancora d'adattare 

//**********CELL 0DS************

  
    bool ExportCell0Ds(const string& outputFilePath,
					   const size_t& n,
					   const unsigned int* const& v1,
					   const unsigned int* const& v2)
		ofstream file;
		file.open(outputFilePath);

		if (file.fail())
		{
			cerr<< "file open failed"<< endl;
			return false;
		}

		file << "# Size of the two vectors"<< endl;
		file << n << endl;

		file << "# vector 1"<< endl;
		for (unsigned int i = 0; i < n; i++)
			file << (i != 0 ? " " : "") << v1[i];
		file << endl;

		file << "# vector 2 "<< endl;
		for (unsigned int i = 0; i < n; i++)
			file << (i != 0 ? " " : "") << v2[i];
		file << endl;

    // Close File
    file.close();

    return true;
		/*
		
        ofstream file("./Cell0Ds.txt");

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

*/
    }

// *********************************** CELL 1 DS ***************************

    bool ExportCell1Ds(PolygonalMesh& mesh)
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

bool ExportCell2Ds(PolyhedrallMesh& mesh){
	
	
	
	
	
	return true;
}

//************************************************* CELL 3DS********************************************

bool ExportCell3Ds(PolyhedrallMesh& mesh){
	
	
	
	return true;
}