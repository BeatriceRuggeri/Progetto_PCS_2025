#include <iostream>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;


// cmake -S ./ -B ./Debug -DCMAKE_BUILD_TYPE="Debug"

int main()
{
    int P;
	int Q;
	int B;
	int C;
	int id_vertice1;
	int id_vertice2;

    cout << "Inserisci un valore per P: ";
    cin >> P;
	
	if (P!=3) {
		cerr << "Valore di P non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per P: ";
		cin >> P;
	}
	else {
		cout << "Hai inserito: " << P << endl;
	}
	
	cout << "Inserisci un valore per Q: ";
    cin >> Q;
	
	if (Q!=3 & Q!= 4 & Q!=5) {
		cerr << "Valore di Q non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per Q: ";
		cin >> Q;
	}
	else {
		cout << "Hai inserito: " << Q << endl;
	
	
	cout << "Inserisci un valore per B: ";
    cin >> B;
	if (B < 1) {
		cerr << "Valore di B non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per B: ";
		cin >> B;
	}
	cout << "Inserisci un valore per C: ";
    cin >> C;
	if (C < 1) {
		cerr << "Valore di C non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per C: ";
		cin >> C;
	}
	
	vector<int> Quadrupla = {P, Q, B, C} 
	cout << Quadrupla << endl;
	
	/*
	cout << "Inserisci l'id del vertice 1: ";
    cin >> id_vertice1;
	cout << "Inserisci l'id del vertice 2: ";
    cin >> id_vertice2;
	
    int max_nodo = massimo_nodo(cubo);

    int distanza = bfs(cubo, id_vertice1, id_vertice2, max_nodo);

    if (distanza != -1)
        cout << "Cammino minimo da " << start << " a " << end << ": " << distanza << " passi\n";
    else
        cout << "Nessun cammino trovato.\n";

    return 0;
	*/
	
	/// Per visualizzare online le mesh:
    /// 1. Convertire i file .inp in file .vtu con https://meshconverter.it/it
    /// 2. Caricare il file .vtu su https://kitware.github.io/glance/app/

    Gedim::UCDUtilities utilities;
    {
        vector<Gedim::UCDProperty<double>> cell0Ds_properties(1);

        cell0Ds_properties[0].Label = "Marker";
        cell0Ds_properties[0].UnitLabel = "-";
        cell0Ds_properties[0].NumComponents = 1;

        vector<double> cell0Ds_marker(mesh.NumCell0Ds, 0.0);
        for(const auto &m : mesh.MarkerCell0Ds)
            for(const unsigned int id: m.second)
                cell0Ds_marker.at(id) = m.first;

        cell0Ds_properties[0].Data = cell0Ds_marker.data();

        utilities.ExportPoints("./Cell0Ds.inp",
                               mesh.Cell0DsCoordinates,
                               cell0Ds_properties);
    }

    {
        vector<Gedim::UCDProperty<double>> cell1Ds_properties(1);

        cell1Ds_properties[0].Label = "Marker";
        cell1Ds_properties[0].UnitLabel = "-";
        cell1Ds_properties[0].NumComponents = 1;

        vector<double> cell1Ds_marker(mesh.NumCell1Ds, 0.0);
        for(const auto &m : mesh.MarkerCell1Ds)
            for(const unsigned int id: m.second)
                cell1Ds_marker.at(id) = m.first;

        cell1Ds_properties[0].Data = cell1Ds_marker.data();

        utilities.ExportSegments("./Cell1Ds.inp",
                                 mesh.Cell0DsCoordinates,
                                 mesh.Cell1DsExtrema,
                                 {},
                                 cell1Ds_properties);
    }
    

    return 0;
}

	
	
	/*
	int Matrix PQ[5][2] = {{3, 3}, {3, 4}, {4, 3}, {5, 3}, {3, 5}};
	
	for (i= 0; i<=5; i++)
	{	
		if (PQ[i]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE1(P, Q, b, c)
		}	
		if (PQ[i][1]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE2 (P, Q, b, c)
		}
	}
	*/
	//trial
	
	PolyhedralLibrary::PolyhedralMesh mesh;
	/*
	bool success_edges_cube = ExportCell1Ds("C_Cell1Ds.txt",mesh.cube_edges,12);
	
	if(!success_edges_cube){
		cerr<<"Export cube's edges failed."<<endl;
		return 1;
	}
	
	bool success_tetrahedron = ExportCell0Ds("T_Cell0Ds.txt",mesh.Vert_tetrahedron,4);
	
	if(!success_tetrahedron){
		cerr<<"Export tetrahedron failed."<<endl;
		return 1;
	}
	
	
	bool success_cube = ExportCell0Ds("C_Cell0Ds.txt",mesh.Vert_cube,8);
	
	if(!success_cube){
		cerr<<"Export cube failed."<<endl;
		return 1;
	}
	
	bool success_octahedron = ExportCell0Ds("O_Cell0Ds.txt",mesh.Vert_octahedron,6);
	
	if(!success_octahedron){
		cerr<<"Export octahedron failed."<<endl;
		return 1;
	}
	
	bool success_dodecahedron = ExportCell0Ds("D_Cell0Ds.txt",mesh.Vert_dodecahedron,20);
	
	if(!success_dodecahedron){
		cerr<<"Export dodecahedron failed."<<endl;
		return 1;
	}
	
	bool success_icosahedron = ExportCell0Ds("I_Cell0Ds.txt",mesh.Vert_icosahedron,12);
	
	if(!success_icosahedron){
		cerr<<"Export icosahedron failed."<<endl;
		return 1;
	}
	*/			
