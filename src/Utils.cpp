#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
using namespace PolyhedralLibrary;



namespace PolyhedralLibrary
{
	/*
    // ***************************************************************************
    bool ExportMesh(PolygonalMesh& mesh)
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
    
    */

//funzione che controlla che i vertici siano sulla sfera di raggio 1 centrata nell'origine. 

/*   
bool ControlloSfera(double x, double y, double z){
	if (x*x +y*y +z*z ==1)
	else 
		cerr << "Vertice non idoneo" << endl;
	return true;
}
*/
    
//*********************************************** Esportazione ********************* : Ancora d'adattare 

//**********CELL 0DS************

  
  bool ExportCell0Ds(const string& outputFilePath, const double vertici[][3],int n)
   {
	
		ofstream file;
		file.open(outputFilePath);

		if (file.fail())
		{
			cerr<< "file open failed"<< endl;
			return false;
		}
		
		file <<"ID "<<"X "<<"Y "<<"Z "<<endl;
		
		for (int i=0; i<n; i++)
		{
			
		   double x = vertici[i][0];
           double y = vertici[i][1];
           double z = vertici[i][2];
	
			
		file << i <<";"<<x<<";"<<y<<";"<<z<<endl;
	    }
		


    // Close File
    file.close();
    
    cout<< "Saved: "<< outputFilePath << endl;

    return true;
   }
/*
// *********************************** CELL 1 DS ***************************

    bool ExportCell1Ds(const string& outputFilePath)
    {
        ofstream file;
		file.open(outputFilePath);

		if (file.fail())
		{
			cerr<< "file open failed"<< endl;
			return false;
		}

		unsigned int n = 7;
		
		file << "id " << "id_start" << "id_end " << endl;
		
		for (i=0; i<=n; i++)
		{
			//condizioni su inzio e file del lato
			double id_start = rand()%n;
			double id_end = rand()%n;
				
			if (id_end-id_start!=0){
				file << i << id_start << id_end << endl;
			}
		}

    // Close File
    file.close();

    return true;
    }
*/

//************************************************* CELL 2DS********************************************

/*
bool ExportCell2Ds(PolyhedrallMesh& mesh){
	
	
	
	
	
	return true;
 }

//************************************************* CELL 3DS********************************************

bool ExportCell3Ds(PolyhedrallMesh& mesh){
	
	
	
	return true;
 }
*/

/*
//************************************** CLASSE1 *****************
bool CLASSE1(...) {
	
	//fa le modifiche che deve avendo (p, q, b, c)
	...
	// fa le esportazioni per ogni poliedro
	
	ExportCell0Ds
	ExportCell1Ds
	ExportCell2Ds
	ExportCell3Ds
	
	return true;
 }

*/

/*
//************************************** CLASSE2 *****************
bool CLASSE2(...) {
	
	//fa le modifiche che deve avendo (p, q, b, c)
	...
	// fa le esportazioni per ogni poliedro
	
	ExportCell0Ds
	ExportCell1Ds
	ExportCell2Ds
	ExportCell3Ds
	
	return true;
 }
*/

}