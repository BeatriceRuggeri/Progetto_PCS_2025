#include "Utils.hpp"          
#include "PolyhedralMesh.hpp" 
#include <iostream>           // Per std::cout, std::cerr
#include <vector>             // Per std::vector
#include <queue>              // Per std::priority_queue
#include <map>                // Per std::map
#include <limits>             // Per std::numeric_limits
#include <cmath>              // Per std::sqrt (usato in calculateDistanceById)
#include <algorithm>          // Per std::min, std::max
#include <tuple>              // Per std::tie o destructuring assignment (C++17)
#include <utility>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen; 

namespace PolyhedralLibrary{

	// IMPLEMENTAZIONE DEL COSTRUTTORE DI SHORTESTPATHRESULT
	minimo_output::minimo_output(unsigned int nEdges, double len)
		: num_lati(nEdges), length(len)
	{}

	// Funzione per calcolare la distanza euclidea tra due punti.
	double euclidian_dist(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2) {
		
		VectorXd A = mesh.Cell0DsCoordinates.col(id1);
		VectorXd B = mesh.Cell0DsCoordinates.col(id2);

		return (A - B).norm(); // Calcola la norma (distanza euclidea)
	}


	// Funzione per calcolare la matrice di adiacenza
	MatrixXi adj_M(const PolyhedralMesh& mesh) {
		const unsigned int Num_vertici = mesh.Cell0DsCoordinates.cols();
		MatrixXi Matrice_adiacenza = MatrixXi::Zero(Num_vertici, Num_vertici);

		// Itera su tutti i lati (Cell1Ds) della mesh	
		for (unsigned int i = 0; i < mesh.Cell1DsId.size(); ++i) {

			unsigned int V1 = mesh.Cell1DsExtrema(i, 0);
			unsigned int V2 = mesh.Cell1DsExtrema(i, 1);
		
			Matrice_adiacenza(V1, V2) = 1;
			Matrice_adiacenza(V2, V1) = 1;
		}

    return Matrice_adiacenza;
	}

	minimo_output minimo_F_Dijkstra(
		PolyhedralMesh& mesh, const MatrixXi& Matrice_adiacenza,
		const unsigned int startVertexId,
		const unsigned int endVertexId) 
		{
			const unsigned int Num_vertici = mesh.Cell0DsCoordinates.cols(); // Numero totale di vertici
			const unsigned int Num_lati_mesh = mesh.Cell1DsId.size();     // Numero totale di lati nella mesh
    
			if (startVertexId >= Num_vertici) {
				cerr << "Il vertice di partenza (" << startVertexId << ") è fuori dal range [0, " << Num_vertici - 1 << "]." << endl;
				// si restituisce un risultato nullo per indicare un errore
				return minimo_output(0, 0.0);
			}
			if (endVertexId >= Num_vertici) {
				cerr << "Il vertice di arrivo (" << endVertexId << ") è fuori dal range [0, " << Num_vertici - 1 << "]." << endl;
				// Restituisci un risultato vuoto per indicare un errore
				return minimo_output(0, 0.0);
			}

			// Inizializza il risultato
			minimo_output risultato(0, 0.0);

			// Esaminiamo il caso banale in cui partenza e arrivo coincidono
			if (startVertexId == endVertexId) {
				cout << "Partenza e arrivo coincidono e ciò comporta un cammino nullo." << endl;
				mesh.Cell0DsMarker.assign(Num_vertici, 0); 
				mesh.Cell1DsMarker.assign(Num_lati_mesh, 0);
				mesh.Cell0DsMarker[startVertexId] = 1; // Marca il vertice sulla mesh
				return risultato;
			}

			// Mappa per collegare una coppia di vertici (tramite INDICI) all'ID del lato e alla sua lunghezza.
			map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> mappa_lati;	
    
			for (unsigned int i = 0; i < Num_lati_mesh; ++i) {
				unsigned int A = mesh.Cell1DsExtrema(i, 0);
				unsigned int B = mesh.Cell1DsExtrema(i, 1);

			double lunghezza = euclidian_dist(mesh, A, B);

			// Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
			mappa_lati[{min(A, B), max(A, B)}] = {mesh.Cell1DsId[i], lunghezza};
			}

			// Variabili Dijkstra
	
			// dist[i] = distanza minima nota dal vertice di partenza a i
			vector<double> dist(Num_vertici, numeric_limits<double>::infinity());
	
			// predVertex[i] = indice del vertice precedente nel cammino minimo a i
			vector<unsigned int> predVertex(Num_vertici, -1); // Usiamo -1 per indicare nessun predecessore
	
			// predEdge[i] = ID del Cell1D usato per raggiungere i dal suo predecessore
			vector<unsigned int> predEdge(Num_vertici, -1);

			// visitato[i] = true se il cammino più breve a 'i' è stato finalizzato
			vector<bool> visitato(Num_vertici, false); 

			// Coda di priorità = memorizza coppie {distanza, indice_vertice}
			// std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
			using QueueElementi = pair<double, unsigned int>;
			priority_queue<QueueElementi, vector<QueueElementi>, greater<>> pq;
	
			// priority_queue<QueueElementi, vector<QueueElementi>, greater<QueueElem>> pq;

			// Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
			dist[startVertexId] = 0.0;
			pq.push({0.0, startVertexId});

			// algoritmo Dijkstra
			while (!pq.empty()) {
			// Estrai il vertice 'u' con la distanza minima corrente dalla coda
				auto [current_distance, u] = pq.top(); // accedo all'elemento
				pq.pop(); // rimuovo l'elemento

			// Se il vertice è già stato visitato, significa che abbiamo già trovato il cammino più breve per esso ==> ignoriamo questa istanza.
				if (visitato[u]) {
					continue;
				}
				visitato[u] = true; // Marca il vertice come visitato/finalizzato

			// Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
				if (u == endVertexId) {
					break;
				}

			// Itera su tutti i possibili vicini 'v' di 'u' usando la matrice di adiacenza
				for (unsigned int v = 0; v < Num_vertici; ++v) {
					// Se c'è un lato tra u e v (adjMatrix(u, v) == 1)
					if (Matrice_adiacenza(u, v) == 1) {
						// Recupera le informazioni sul lato (ID reale e lunghezza/peso)
						auto it_lato = mappa_lati.find({min(u, v), max(u, v)});
						if (it_lato == mappa_lati.end()) {
							cerr << "C'è un lato tra gli indici (" << u << "," << v << ") presente nella matrice di adiacenza ma non trovato nella mappa dei lati.\n";
							continue;
						}
					double peso = it_lato->second.second; // Peso del lato (lunghezza euclidea)

					// Operazione di "rilassamento":
					// Se la distanza calcolata a 'v' passando per 'u' è minore della distanza attuale di 'v'
					if (dist[u] + peso < dist[v]) {
						dist[v] = dist[u] + peso; // Aggiorna la distanza minima per 'v'
						predVertex[v] = u;          // Imposta 'u' come predecessore di 'v'
						predEdge[v] = it_lato->second.first; // Memorizza l'ID reale del lato usato
						pq.push({dist[v], v});      // Inserisci 'v' nella coda di priorità con la nuova distanza
					}
				}
			}
		}
	
		// Se la distanza al vertice di arrivo è ancora infinito, significa che non è stato trovato alcun cammino
		if (dist[endVertexId] == numeric_limits<double>::infinity()) {
			cout << "Non è presente nessun cammino tra il vertice " << startVertexId
                  << " e il vertice " << endVertexId << endl;
			mesh.Cell0DsMarker.assign(Num_vertici, 0);
			mesh.Cell1DsMarker.assign(Num_lati_mesh, 0);
			return risultato;
		}
	
	
		// Ricostruzione del cammino e calcolo delle statistiche
		// Ricostruisci il cammino a ritroso dal vertice di arrivo al vertice di partenza
	
		// Inizializzo i marker della mesh a 0
		mesh.Cell0DsMarker.assign(Num_vertici, 0); // Inizializza i marker a 0
		mesh.Cell1DsMarker.assign(Num_lati_mesh, 0);
	
		risultato.length = dist[endVertexId]; // distanza totale

		unsigned int id_vert_preso = endVertexId; // Partiamo dall'indice del vertice di arrivo
	
		while (id_vert_preso != startVertexId) {
			// Marca il vertice corrente come parte del cammino minimo
			mesh.Cell0DsMarker[id_vert_preso] = 1;

			unsigned int id_vert_prec = predVertex[id_vert_preso]; // Indice del vertice precedente nel cammino
			unsigned int id_lato_preso = predEdge[id_vert_preso];      // ID del lato usato per raggiungere current_idx

			mesh.Cell1DsMarker[id_lato_preso] = 1;
			risultato.num_lati++;
        
			id_vert_preso = id_vert_prec; // Spostati al vertice precedente e continua la ricostruzione
		}
	
		// Marca anche il vertice di partenza come parte del cammino
		mesh.Cell0DsMarker[startVertexId] = 1;	
	
		return risultato;
	}

}