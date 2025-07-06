#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>

#include <limits> 
#include <Eigen/Dense>
#include <set> 
#include <cmath>

namespace PolyhedralLibrary {

	// funzione che assegna una flag (come fosse un marker) ai vertici che devono rimanere (non duplicati)
	void v_duplicates_rm(PolyhedralMesh& meshTriangolata)
	{
		double tol = 1e-12;
		unsigned int Max_Flag = numeric_limits<unsigned int>::max();
		size_t n = meshTriangolata.Cell0DsCoordinates.cols(); 

		// inizializzo un vettore per tenere traccia del reindirizzamento finale degli ID dei vertici.
		// imposto che all'inizio ogni vertice punta a se stesso.
		vector<unsigned int> id_reindirizzati(n);
		for (unsigned int i = 0; i < n; ++i) {
			id_reindirizzati[i] = i;
		}
	
		// inizializzo un vettore che tiene traccia di quali vertici sono potenziali non duplicati. 
		// impongo che se un vertice viene identificato come duplicato di un altro, la sua flag viene impostata a false
		vector<bool> miglior_candidato(n, true);

		for (size_t i = 0; i < n; ++i) {
			// Se il vertice 'i' è già stato marcato per essere tenuto ==> Max_Flag=TRUE
			// Se è stato reindirizzato da un vertice precedente nel ciclo significa che è un duplicato ==> salta.
			if (meshTriangolata.Cell0DsFlag[i][0] == Max_Flag || !miglior_candidato[i])
				continue;

			for (size_t j = i + 1; j < n; ++j) {
				if (meshTriangolata.Cell0DsFlag[j][0] == Max_Flag || !miglior_candidato[j])
					continue;

				bool FattorComune = false;
				for (unsigned int Flag_i : meshTriangolata.Cell0DsFlag[i]) {
					for (unsigned int Flag_j : meshTriangolata.Cell0DsFlag[j]) {
						if (Flag_i == Flag_j) {
							FattorComune = true;
							break;
						}
					}
					if (FattorComune) break;
				}

				if (FattorComune) {
					if ((meshTriangolata.Cell0DsCoordinates.col(i) - meshTriangolata.Cell0DsCoordinates.col(j)).norm() < tol) {
						// Il vertice 'i' è un duplicato di 'j', vogliamo mantenere 'j' e reindirizzare 'i' a 'j'.
						// Marchiamo 'i' come duplicato
						miglior_candidato[i] = false;
						// Reindirizziamo 'i' a 'j'.
						// Ora tutti i vertici che precedentemente puntavano a 'i' devono puntare a 'j' e 'i' stesso punterà a 'j'.
                    
						unsigned int ind_i = i; //indirizzo di i
						while (id_reindirizzati[ind_i] != ind_i) {
							ind_i = id_reindirizzati[ind_i];
						}
						unsigned int ind_j = j; //indirizzo di j
						while (id_reindirizzati[ind_j] != ind_j) {
							ind_j = id_reindirizzati[ind_j];
						}

						if (ind_i != ind_j){ 
							id_reindirizzati[ind_i] = ind_j;
						}
						// vertice 'i' (e i suoi precedenti duplicati) puntano a 'j' o alla sua flag
						id_reindirizzati[i] = ind_j;

						// Assegnamo le coordinate di 'j' a 'i'
						meshTriangolata.Cell0DsCoordinates.col(i) = meshTriangolata.Cell0DsCoordinates.col(j);
					}
				}
			}
		}
	    
		// Applichiamo il reindirizzamento finale a tutti i vertici e dobbiamo propagare le catene di reindirizzamento.
		for (unsigned int i = 0; i < n; ++i) {
			unsigned int id_attuale = i;
			while (id_reindirizzati[id_attuale] != id_attuale) {
				id_attuale = id_reindirizzati[id_attuale];
			}
			id_reindirizzati[i] = id_attuale; // Imposta l'ultimo reindirizzamento per i
		}

		// Aggiorniamo le strutture create all'inizio in base al reindirizzamento finale
		for (int edgeId = 0; edgeId < meshTriangolata.Cell1DsExtrema.rows(); ++edgeId) {
			meshTriangolata.Cell1DsExtrema(edgeId, 0) = id_reindirizzati[meshTriangolata.Cell1DsExtrema(edgeId, 0)];
			meshTriangolata.Cell1DsExtrema(edgeId, 1) = id_reindirizzati[meshTriangolata.Cell1DsExtrema(edgeId, 1)];
		}
    
		for (unsigned int faceId = 0; faceId < meshTriangolata.Cell2DsId.size(); ++faceId) {
			for (unsigned int& vertexOriginalId : meshTriangolata.Cell2DsVertices[faceId]) {
				vertexOriginalId = id_reindirizzati[vertexOriginalId];
			}
		}

		// Aggiorniamo Cell0DsId e Cell0DsFlag in base al reindirizzamento.
		// I vertici che hanno come Flag miglior_candidato, avranno il loro id_reindirizzati[i] == i mentre quelli duplicati avranno id_reindirizzati[i] != i.
		for (unsigned int i = 0; i < n; ++i) {
			meshTriangolata.Cell0DsId[i] = id_reindirizzati[i];
			if (id_reindirizzati[i] == i) {
				// Questo vertice è un miglior_candidato o non è stato reindirizzato
				meshTriangolata.Cell0DsFlag[i] = {Max_Flag};
			} 
		}
	}

	// funzione che assegna una flag (come fosse un marker) ai lati/spigoli che devono rimanere (non duplicati)
		void l_duplicates_rm(PolyhedralMesh& meshTriangolata)
	{
		double tol = 1e-12;
		unsigned int Max_Flag = numeric_limits<unsigned int>::max();
		size_t n = meshTriangolata.Cell1DsExtrema.rows();
    
		if (meshTriangolata.Cell1DsOriginalFlag.empty() || meshTriangolata.Cell1DsOriginalFlag.size() != n) {
			meshTriangolata.Cell1DsOriginalFlag.resize(n);
		}
    
		for (size_t i = 0; i < n; ++i) {
			meshTriangolata.Cell1DsOriginalFlag[i] = (meshTriangolata.Cell1DsFlag[i] == Max_Flag);
		}

		// inizializzo un vettore per tenere traccia del reindirizzamento finale degli ID dei lati, impostando che ogni lato punta a se stesso.
		vector<unsigned int> id_reindirizzati(n);
		for (unsigned int i = 0; i < n; ++i) {
			id_reindirizzati[i] = i;
		}

		// inizializzo un vettore che tiene traccia di quali lati sono potenziali non duplicati. 
		// Se un lato viene identificato come duplicato di un altro ==> la sua flag miglior candidato viene impostata a false
		vector<bool> miglior_candidato(n, true);

		for (size_t i = 0; i < n; ++i) {
			// Se il lato corrente 'i' è già stato marcato per essere tenuto ==> maxFlag= TRUE
			// se è stato reindirizzato da un lato precedente nel ciclo ==> salta.
			if (meshTriangolata.Cell1DsFlag[i] == Max_Flag || !miglior_candidato[i])
				continue;

			for (size_t j = i + 1; j < n; ++j) {
				if (meshTriangolata.Cell1DsFlag[j] == Max_Flag || !miglior_candidato[j])
					continue;

				if (meshTriangolata.Cell1DsFlag[i] == meshTriangolata.Cell1DsFlag[j]) {
					int i0 = meshTriangolata.Cell1DsExtrema(i, 0); // Primo estremo del lato i
					int i1 = meshTriangolata.Cell1DsExtrema(i, 1); // Secondo estremo del lato i
					int j0 = meshTriangolata.Cell1DsExtrema(j, 0); // Primo estremo del lato j
					int j1 = meshTriangolata.Cell1DsExtrema(j, 1); // Secondo estremo del lato j

					// Controlliamo se i vertici estremi corrispondono in ordine diretto o inverso
					bool corrispondenza_diretta = ((meshTriangolata.Cell0DsCoordinates.col(i0) - meshTriangolata.Cell0DsCoordinates.col(j0)).norm() < tol &&
												(meshTriangolata.Cell0DsCoordinates.col(i1) - meshTriangolata.Cell0DsCoordinates.col(j1)).norm() < tol);

					bool corrispondenza_inversa = ((meshTriangolata.Cell0DsCoordinates.col(i0) - meshTriangolata.Cell0DsCoordinates.col(j1)).norm() < tol &&
												(meshTriangolata.Cell0DsCoordinates.col(i1) - meshTriangolata.Cell0DsCoordinates.col(j0)).norm() < tol);

					if (corrispondenza_diretta || corrispondenza_inversa) {
						//  Il lato 'i' è un duplicato del lato 'j' ==> manteniamo 'j' e reindirizziamo 'i' a 'j'.

						miglior_candidato[i] = false; // identifichiamo 'i' come duplicato	.
	
						unsigned int ind_i = i; //indirizzo di i
						while (id_reindirizzati[ind_i] != ind_i) {
							ind_i = id_reindirizzati[ind_i];
						}
						unsigned int ind_j = j; //indirizzo di j
						while (id_reindirizzati[ind_j] != ind_j) {
							ind_j = id_reindirizzati[ind_j];
						}

						if (ind_i != ind_j) {
							id_reindirizzati[ind_i] = ind_j;
						}
						// imponiamo che il lato 'i' stesso (e i suoi precedenti duplicati) puntino a 'j' o alla sua flag.
						id_reindirizzati[i] = ind_j; // Imposta il reindirizzamento diretto per 'i'

						//meshTriangulated.Cell1DsExtrema.row(i) = meshTriangulated.Cell1DsExtrema.row(j);
					}
				}
			}
		}

		for (unsigned int i = 0; i < n; ++i) {
			unsigned int id_attuale = i;
			while (id_reindirizzati[id_attuale] != id_attuale) {
				id_attuale = id_reindirizzati[id_attuale];
				id_reindirizzati[i] = id_attuale;
			}
			id_reindirizzati[i] = id_attuale;
		}

		// Aggiorniamo le strutture definite in base al reindirizzamento, aggiornando i riferimenti ai lati nelle facce (Cell2DsEdges)
		for (unsigned int faceId = 0; faceId < meshTriangolata.Cell2DsId.size(); ++faceId) {
			for (unsigned int& edgeOriginalId : meshTriangolata.Cell2DsEdges[faceId]) {
				edgeOriginalId = id_reindirizzati[edgeOriginalId];
			}
		}

		// Aggiorniamo Cell1DsId e Cell1DsFlag con i risultati del reindirizzamento.
		for (unsigned int i = 0; i < n; ++i) {
			meshTriangolata.Cell1DsId[i] = id_reindirizzati[i];
			if (id_reindirizzati[i] == i) {
				meshTriangolata.Cell1DsFlag[i] = Max_Flag;
			} else {
				// Questo lato è un duplicato e punterà a un master
			}
		}
	}

	void meshClean(const PolyhedralMesh& meshTriangolata, PolyhedralMesh& meshOutput, const vector<int>& dimension)
	{
		unsigned int Max_Flag = numeric_limits<unsigned int>::max();
		meshOutput.Cell0DsId.resize(dimension[0]);
		meshOutput.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
		
		meshOutput.Cell1DsId.resize(dimension[1]);
		meshOutput.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
		meshOutput.Cell1DsOriginalFlag.resize(dimension[1]);
		
		meshOutput.Cell2DsId.resize(dimension[2]);
		meshOutput.Cell2DsVertices.resize(dimension[2]);
		meshOutput.Cell2DsEdges.resize(dimension[2]);
		
		//------------------------------------------------------------------------------------------------------------------
		// VERTICI
    	map<unsigned int, unsigned int> mappa_conv_vert; // Mappa per convertire i vecchi ID dei vertici con i nuovi ID dei vertici
    	vector<unsigned int> vertici_originali; // Tiene conto degli id originali dei vertici unici
    	
    	unsigned int A = 0; // Contatore per i nuovi ID dei vertici
    	for (unsigned int i = 0; i < meshTriangolata.Cell0DsCoordinates.cols(); ++i) {
			if (meshTriangolata.Cell0DsFlag[i][0] == Max_Flag) { 
				mappa_conv_vert[i] = A; // Mappa il vecchio ID unico al nuovo
				vertici_originali.push_back(i); // Salva il vecchio ID per recuperare le coordinate dopo
				A++; // Incrementa il contatore del nuovo ID
			}
		}

		// Compiliamo le strutture finali per i vertici
		for (unsigned int i = 0; i < A; ++i) {
			meshOutput.Cell0DsId[i] = i; // Nuovi ID consecutivi
			meshOutput.Cell0DsCoordinates.col(i) = meshTriangolata.Cell0DsCoordinates.col(vertici_originali[i]);
		}
		
		//------------------------------------------------------------------------------------------------------------------------
		// SPIGOLI
		map<unsigned int, unsigned int> mappa_conv_lati; // Mappa per convertire i vecchi ID degli spigoli con i nuovi ID degli spigol
		vector<pair<unsigned int, unsigned int>> estremi_lati_tmp; // Vettore temporaneo per gli estremi degli spigoli (con i nuovi ID vertici)
		vector<bool> lati_originali_flag_tmp; // Temp per il nuovo flag
		
		unsigned int B = 0; // Contatore per i nuovi ID degli spigoli
		for (unsigned int i = 0; i < meshTriangolata.Cell1DsExtrema.rows(); ++i) {
			if (meshTriangolata.Cell1DsFlag[i] == Max_Flag) {
				unsigned int id_v1_vecchio = meshTriangolata.Cell1DsExtrema(i, 0);
				unsigned int id_v2_vecchio = meshTriangolata.Cell1DsExtrema(i, 1);
	
				// Ottieniamo i nuovi ID dei vertici.
				// Se uno dei vertici originali non è stato mantenuto ==> lo spigolo non è valido.
				if (mappa_conv_vert.count(id_v1_vecchio) > 0 && mappa_conv_vert.count(id_v2_vecchio) > 0) { //restituisce 1 se la chiave esiste, 0 altrimenti.
					unsigned int id_v1_nuovo = mappa_conv_vert[id_v1_vecchio];
					unsigned int id_v2_nuovo = mappa_conv_vert[id_v2_vecchio];
	
					estremi_lati_tmp.push_back({min(id_v1_nuovo, id_v2_nuovo), max(id_v1_nuovo, id_v2_nuovo)});
					mappa_conv_lati[i] = B; // Mappa il vecchio ID spigolo al nuovo
					lati_originali_flag_tmp.push_back(meshTriangolata.Cell1DsOriginalFlag[i]);
					B++;
				} 
			}
		}
		
		// Compila le strutture finali per gli spigoli
		for (unsigned int i = 0; i < B; ++i) {
			meshOutput.Cell1DsId[i] = i; // Nuovi ID consecutivi
			meshOutput.Cell1DsExtrema(i, 0) = estremi_lati_tmp[i].first;
			meshOutput.Cell1DsExtrema(i, 1) = estremi_lati_tmp[i].second;
			
			meshOutput.Cell1DsOriginalFlag[i] = lati_originali_flag_tmp[i];
		}
		
		//----------------------------------------------------------------------------------------------------------------------------------
		// FACCE
		vector<vector<unsigned int>> vert_facce_tmp; // Vettore temporaneo per i vertici delle facce (con i nuovi ID vertici)
		vector<vector<unsigned int>> spigoli_facce_tmp; // Vettore temporaneo per gli spigoli delle facce (con i nuovi ID spigoli)
		
		unsigned int C = 0; // Contatore per i nuovi ID delle facce
		for (unsigned int i = 0; i < meshTriangolata.Cell2DsId.size(); ++i) { // Copia e aggiorna i vertici della faccia
			vector<unsigned int> nuovi_vert; //indichiamo i nuovi vertici della faccia corrente
			bool vertici_validi = true; //indichiamo come variabile booleana per indicare se i vertici nuovi sono validi per le facce o meno
			for (unsigned int vertici_vecchi : meshTriangolata.Cell2DsVertices[i]) {
				if (mappa_conv_vert.count(vertici_vecchi) > 0) { // Se il vertice è stato mantenuto
					nuovi_vert.push_back(mappa_conv_vert[vertici_vecchi]);
				} else {
					vertici_validi = false; // Vertice non valido, la faccia è invalida
					break;
				}
			}
	
			// Copia e aggiorna gli spigoli della faccia
			vector<unsigned int> nuovi_lati; //vettore che va a registrare tutti i nuovi lati della faccia corrente
			bool lati_validi = true; //indichiamo come variabile booleana per indicare se i lati nuovi sono validi per le facce o meno
			for (unsigned int lati_vecchi : meshTriangolata.Cell2DsEdges[i]) {
				if (mappa_conv_lati.count(lati_vecchi) > 0) { // Se lo spigolo è stato mantenuto
					nuovi_lati.push_back(mappa_conv_lati[lati_vecchi]);
				} else {
					lati_validi = false; // Spigolo non valido, la faccia è invalida
					break;
				}
			}
			
			// Se la faccia è valida (tutti i suoi vertici e spigoli sono stati mantenuti e riassegnati)
			if (vertici_validi && lati_validi) {
				vert_facce_tmp.push_back(nuovi_vert);
				spigoli_facce_tmp.push_back(nuovi_lati);
				C++;
			}
		}
	
		// compila le strutture finali per le facce
		for (unsigned int i = 0; i < C; ++i) {
			meshOutput.Cell2DsId[i] = i; // Nuovi ID consecutivi
			meshOutput.Cell2DsVertices[i] = vert_facce_tmp[i];
			meshOutput.Cell2DsEdges[i] = spigoli_facce_tmp[i];
		}
		
		
	}
}