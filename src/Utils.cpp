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
// *********************************** CELL 1 DS ***************************

    bool ExportCell1Ds(const string& outputFilePath, const vector<Edge>& edges)
    {   
        ofstream file;
		file.open(outputFilePath);

		if (file.fail())
		{
			cerr<< "file open failed"<< endl;
			return false;
		}
		
		file << "Id " << "Id_start " << "Id_end " << endl;
		
		for (const auto& edge : edges)
		{
			file << edge.id << ";" <<  edge.origin << ";" << edges.end << endl;
		}

    // Close File
    file.close();

    return true;
    }


/*
int main() {
    // Lista di vertici (ognuno è un punto 2D)
    std::vector<std::vector<float>> vertices = {
        {0.0f, 0.0f},   // Vertice 0
        {1.0f, 0.0f},   // Vertice 1
        {1.0f, 1.0f},   // Vertice 2
        {0.0f, 1.0f}    // Vertice 3
    };

    // Vettore per memorizzare gli edge
    std::vector<Edge> edges;

    // Definiamo gli edge (i lati del quadrato)
    int numVertices = vertices.size();
    for (int i = 0; i < numVertices; ++i) {
        // Ogni edge collega il vertice i con il successivo, e l'ultimo con il primo
        int startIdx = i;
        int endIdx = (i + 1) % numVertices;  // Per creare il ciclo, % numVertices
        edges.push_back(Edge(startIdx, endIdx));
    }

    // Stampa gli edge
    for (const Edge& e : edges) {
        std::cout << "Edge from vertex " << e.start << " to vertex " << e.end << std::endl;
    }

    return 0;
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