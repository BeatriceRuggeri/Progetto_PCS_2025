#include <iostream>
#include <cstdlib> 
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main(int argc, char *argv[]) {

    // Definizione delle variabili per il cammino minimo, inizializzate a valori non validi
    unsigned int id_vertice_inizio = 0;
    unsigned int id_vertice_fine = 0;
    bool calcolatoPercorso = false; // Flag per indicare se calcolare il cammino minimo

    if (argc == 5) { // Caso: solo p, q, b, c
        cout << "in esecuzione: creazione mesh.\n";
        // Nessun cammino minimo da calcolare
        calcolatoPercorso = false;

    } else if (argc == 7) { // Caso: p, q, b, c, id_vertice_inizio, id_vertice_fine
        cout << "in esecuzione: creazione mesh e calcolo cammino minimo.\n";
        id_vertice_inizio = stoul(argv[5]);
        id_vertice_fine = stoul(argv[6]);
        calcolatoPercorso = true;

    } else { // Errore: numero di argomenti non valido
        cerr << "Uso:\n";
        cerr << "  " << argv[0] << " p q b c\n";
        cerr << "  " << argv[0] << " p q b c id_vertice_inizio id_vertice_fine\n";
        return 1;
    }
    
    int p = stoi(argv[1]);
    int q = stoi(argv[2]);
    int b = stoi(argv[3]);
    int c = stoi(argv[4]);

	if (!((p == 3 || q == 3) && (p >= 3 && p <= 5) && (q >= 3 && q <= 5))) {
		cerr << "Attenzione, almeno uno tra p e q deve valere 3 e i valori ammessi per entrambi sono 3, 4 e 5\n";
		return 1;
	}
	
	if (!((b == c) || (b == 0) || (c == 0))) {
		cerr << "Attenzione, si richiede c = 0, oppure b = c.\n";
		return 1;
	}
	
	if ((b == c) && (p!=3)){
		cerr << "Attenzione, richiestro poliedro geodetico per la triangolazione di classe II.\n";
		return 1;
	}
	
    cout << "quaterna digitata: p = " << p << ", q = " << q << ", b = " << b << ", c = " << c << ".\n";
    if (calcolatoPercorso) {
        cout << "Cammino minimo da vertice " << id_vertice_inizio << " a vertice " << id_vertice_fine << ".\n";
    }

	PolyhedralMesh mesh;
	PolyhedralMesh meshOutput;
    
    if (b != c){
	    
		if (p == 3 && q == 3) {
			tetraedro_gen(mesh);
			triangolazione_I(q, b, c, mesh, meshOutput);
			
		} else if (p == 3 && q != 3){
			if (q == 4){
				octaedro_gen(mesh);
			} else {
				icosaedro_gen(mesh);
			}
			triangolazione_I(q, b, c, mesh, meshOutput);
	
		} else if (q == 3 && p!= 3) {
			if (p == 4){
				octaedro_gen(mesh);
			} else {
				icosaedro_gen(mesh);
			}
			swap_val(p, q);
			triangolazione_dual(q, b, c, mesh, meshOutput);
		}
	} else {
		if (p == 3 && q == 3) {
			tetraedro_gen(mesh);
			triangolazione_II(q, b, mesh, meshOutput);
			
		} else if (p == 3 && q != 3){
			if (q == 4){
				octaedro_gen(mesh);
			} else {
				icosaedro_gen(mesh);
			}
			triangolazione_II(q, b, mesh, meshOutput);
		}
	} 

	if (calcolatoPercorso) {
	MatrixXi adjMatrix = adj_M(meshOutput);
        
        minimo_output pathResult = minimo_F_Dijkstra(
            meshOutput,
            adjMatrix,
            id_vertice_inizio,
            id_vertice_fine
        );

        // Stampa i risultati
        if (pathResult.num_lati > 0 || id_vertice_inizio == id_vertice_fine) {
            cout << "\n|||| Cammino Minimo ||||\n";
            cout << "Numero di lati nel cammino: " <<pathResult.num_lati<< endl;
            cout << "Lunghezza totale del cammino: "<<pathResult.length<< endl;

            // Esportazione Paraview con cammino
			proj_sphere(meshOutput);
			ExportParaview(meshOutput);
		
        } else {
            cout << "\nNessun cammino trovato tra il vertice " << id_vertice_inizio
                 << " e il vertice " << id_vertice_fine << ".\n";
        }
	} else {
		// Esportazione Paraview senza cammino
		proj_sphere(meshOutput);
		ExportParaview(meshOutput);
	}

	// Scrittura su TXT
	WriteCell0Ds(meshOutput);
	WriteCell1Ds(meshOutput);
	WriteCell2Ds(meshOutput);
	WriteCell3Ds(meshOutput);

    return 0;
    
}