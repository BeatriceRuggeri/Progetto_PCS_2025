#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>
#include <Eigen/Dense> 
#include <cmath>

using namespace std;
using namespace Eigen;
namespace PolyhedralLibrary {
	
// Funzione per ottenere il baricentro di una faccia
Vector3d baricenter_Fetch(const PolyhedralMesh& mesh_triangolata, const unsigned int id_faccia) {
    Vector3d baricentro = Vector3d::Zero();
    const auto& vertici_faccia = mesh_triangolata.Cell2DsVertices[id_faccia];
    for (unsigned int v_id : vertici_faccia) {
        baricentro += mesh_triangolata.Cell0DsCoordinates.col(v_id);
    }
    return baricentro /= (3.0);
}

// Mappa Spigolo Originale -> Facce che lo Contengono
map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildMap_LF(const PolyhedralMesh& mesh_triangolata) {
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> lato_faccia;

    for (unsigned int id_faccia = 0; id_faccia < mesh_triangolata.Cell2DsId.size(); ++id_faccia) {
        const vector<unsigned int>& faceEdges = mesh_triangolata.Cell2DsEdges[id_faccia]; // spigoli della faccia corrente
        
        for (unsigned int id_originale_lato : faceEdges) {
			unsigned int v1_id = mesh_triangolata.Cell1DsExtrema(id_originale_lato, 0);
			unsigned int v2_id = mesh_triangolata.Cell1DsExtrema(id_originale_lato, 1);
			pair<unsigned int, unsigned int> sortedEdgeVertices = {min(v1_id, v2_id), max(v1_id, v2_id)};
			// Ordiniamo i vertici dello spigolo per avere una chiave univoca nella mappa
			lato_faccia[sortedEdgeVertices].push_back(id_faccia);
        }
    }
    return lato_faccia;
}

// Mappa Vertice Originale -> Facce che lo Contengono
map<unsigned int, vector<unsigned int>> buildMap_VF(const PolyhedralMesh& mesh_triangolata) {
    map<unsigned int, vector<unsigned int>> vertice_faccia;
    for (unsigned int id_faccia = 0; id_faccia < mesh_triangolata.Cell2DsId.size(); ++id_faccia) {
        for (unsigned int id_verticiOriginali : mesh_triangolata.Cell2DsVertices[id_faccia]) {
            vertice_faccia[id_verticiOriginali].push_back(id_faccia);
        }
    }
    return vertice_faccia;
}

// Mappa Vertice Originale -> Spigoli che vi incidono
map<unsigned int, vector<unsigned int>> buildMap_VE(const PolyhedralMesh& mesh_triangolata) {
    map<unsigned int, vector<unsigned int>> vertice_lato;
    for (unsigned int edgeId = 0; edgeId < mesh_triangolata.Cell1DsExtrema.rows(); ++edgeId) {
        unsigned int v1 = mesh_triangolata.Cell1DsExtrema(edgeId, 0);
        unsigned int v2 = mesh_triangolata.Cell1DsExtrema(edgeId, 1);
        vertice_lato[v1].push_back(edgeId);
        vertice_lato[v2].push_back(edgeId);
    }
    return vertice_lato;
}


void dual_calc(PolyhedralMesh& mesh_triangolata, PolyhedralMesh& mesh_duale, map<pair<unsigned int, unsigned int>, vector<unsigned int>> lato_facciaMap)
{
	// VERTICI
	// i vertici del duale sono i baricentri delle facce originali
	mesh_duale.Cell0DsId.resize(mesh_triangolata.Cell2DsId.size());
	mesh_duale.Cell0DsCoordinates = MatrixXd::Zero(3, mesh_triangolata.Cell2DsId.size());
	
	for (unsigned int id_faccia = 0; id_faccia < mesh_triangolata.Cell2DsId.size() ; ++id_faccia){
		mesh_duale.Cell0DsCoordinates.col(id_faccia) = baricenter_Fetch(mesh_triangolata, id_faccia);
        mesh_duale.Cell0DsId[id_faccia] = id_faccia; // L'ID del vertice duale è l'ID della faccia originale
    }
    
    // SPIGOLI
    
    vector<pair<unsigned int, unsigned int>> estremi_lati_duale;
    // Ogni spigolo interno del poliedro originale (condiviso da due facce) genera uno spigolo nel poliedro duale che connette i baricentri di quelle due facce.
    // Useremo un vettore temporaneo e poi lo convertiremo in Eigen::MatrixXi.
        
    for (const auto& entry : lato_facciaMap) {
        const vector<unsigned int>& lato_condiviso_facce = entry.second; // facce
        if (lato_condiviso_facce.size() == 2) {
			unsigned int id_faccia1 = lato_condiviso_facce[0]; // faccia condivisa 1
			unsigned int id_faccia2 = lato_condiviso_facce[1]; // faccia condivisa 2
			pair<unsigned int, unsigned int> dualEdge = {min(id_faccia1, id_faccia2), max(id_faccia1, id_faccia2)};
			// Gli ID delle facce originali diventano gli ID dei vertici del duale
			estremi_lati_duale.push_back(dualEdge);
		}
	}

	// Questa struttura intermedia mi serve perchè non so a prescindere le dimensioni di mesh_duale.Cell1DsExtrema
	// il numero di spigoli duali è pari al numero di spigoli interni della mesh originale
		
	mesh_duale.Cell1DsId.resize(estremi_lati_duale.size());
	mesh_duale.Cell1DsExtrema = MatrixXi::Zero(estremi_lati_duale.size(), 2);
	
	for (unsigned int i = 0; i < estremi_lati_duale.size(); ++i) {
		mesh_duale.Cell1DsId[i] = i;
		mesh_duale.Cell1DsExtrema(i, 0) = estremi_lati_duale[i].first;
		mesh_duale.Cell1DsExtrema(i, 1) = estremi_lati_duale[i].second;
    }
    
    // FACCE
    
    map<pair<unsigned int, unsigned int>, unsigned int> dualEdgeToIdMap;
    // mappa che per ogni coppia di id di facce originali mi associa l'id dello spigolo duale
    for (unsigned int i = 0; i < mesh_duale.Cell1DsId.size(); ++i) {
        unsigned int v1 = mesh_duale.Cell1DsExtrema(i, 0);
        unsigned int v2 = mesh_duale.Cell1DsExtrema(i, 1);
        dualEdgeToIdMap[{min(v1, v2), max(v1, v2)}] = i;
    }
    	 
    map<unsigned int, vector<unsigned int>> vertice_facciaMap = buildMap_VF(mesh_triangolata); 
    map<unsigned int, vector<unsigned int>> vertice_latoMap = buildMap_VE(mesh_triangolata);

    mesh_duale.Cell2DsVertices.resize(vertice_facciaMap.size());
    mesh_duale.Cell2DsId.resize(vertice_facciaMap.size());
    mesh_duale.Cell2DsEdges.resize(vertice_facciaMap.size());

    unsigned int counter_idFacceDuale = 0; // contatore facce duali
    for (const auto& entry : vertice_facciaMap) {
        unsigned int id_verticiOriginali = entry.first; // Il vertice originale che "genera" questa faccia duale
        vector<unsigned int> facceIncidenti = entry.second;  // Le facce originali incidenti a questo vertice
        vector<unsigned int> verticiIncidenti = vertice_latoMap[id_verticiOriginali]; // Gli spigoli originali incidenti a questo vertice

        // Un vertice deve essere incidente ad almeno 3 facce/spigoli per formare una faccia duale chiusa.
        if (facceIncidenti.size() < 3 || verticiIncidenti.size() < 3) {
            cerr << "Attenzione, vertice originale " << id_verticiOriginali << " incidente a meno di 3 facce o spigoli (dunque non può formare una faccia del duale chiusa)" << endl;
            continue; // Salta la creazione della faccia duale
        }

        vector<unsigned int> verticiOrdinati_facciaDuale; // Vertici della faccia duale, ordinati
        vector<unsigned int> latiOrdinati_facciaDuale;           // Spigoli della faccia duale, ordinati

        // Iniziamo il percorso topologico:
        // Scegliamo uno spigolo incidente qualsiasi e una delle facce che lo contiene come punto di partenza.
        unsigned int currentid_originale_lato = verticiIncidenti[0]; // lato di partenza

        // Troviamo le due facce che condividono questo spigolo originale.
        pair<unsigned int, unsigned int> currentEdgeKey = {
            min(mesh_triangolata.Cell1DsExtrema(currentid_originale_lato, 0), mesh_triangolata.Cell1DsExtrema(currentid_originale_lato, 1)),
            max(mesh_triangolata.Cell1DsExtrema(currentid_originale_lato, 0), mesh_triangolata.Cell1DsExtrema(currentid_originale_lato, 1))
        };
        vector<unsigned int> facce_condividenti_lato = lato_facciaMap.at(currentEdgeKey); 
        // Se la currentEdgeKey non esiste nella lato_facciaMap, at() lancerà un'eccezione di tipo std::out_of_range
        // L'operatore [], se la chiave non esiste, inserirà automaticamente una nuova coppia chiave-valore nella mappa, 
        // con il valore predefinito per il tipo di valore

        // Ci assicuriamo che lo spigolo sia condiviso da esattamente due facce
        if (facce_condividenti_lato.size() != 2) {
             cerr << "Errore: Spigolo originale " << currentid_originale_lato << " non condiviso da esattamente due facce" << endl;
             continue; // Salta, non possiamo costruire la faccia duale correttamente
        }

        // Scelgiamo una delle due facce come punto di partenza per il ciclo
        unsigned int startid_faccia = facce_condividenti_lato[0];
        unsigned int idVerticeDuale_precedente = startid_faccia; // Il primo vertice della faccia duale (id della faccia originale)

        // Aggiungiamo il primo vertice duale della faccia che stiamo costruendo
        verticiOrdinati_facciaDuale.push_back(idVerticeDuale_precedente);

        // Percorriamo il ciclo attorno al vertice originale
        for (size_t i = 0; i < verticiIncidenti.size(); ++i) { // Il numero di spigoli incidenti è il numero di lati della faccia duale
            // Troviamo le due facce che condividono l'attuale spigolo `currentid_originale_lato`
            vector<unsigned int> currentEdgeFaces = lato_facciaMap.at(currentEdgeKey);

            unsigned int face1 = currentEdgeFaces[0];
            unsigned int face2 = currentEdgeFaces[1];

            // La faccia "precedente" è idVerticeDuale_precedente. Troviamo la faccia "successiva".
            unsigned int nuovo_id_verticeDuale;
            if (face1 == idVerticeDuale_precedente) {
                nuovo_id_verticeDuale = face2;
            } else if (face2 == idVerticeDuale_precedente) {
                nuovo_id_verticeDuale = face1;
            } else {
                // Questa situazione non dovrebbe accadere se la logica di navigazione è corretta
                cerr << "Errore logico: Faccia corrente non trovata tra le facce che condividono lo spigolo." << endl;
                break;
            }
            
            // Troviamo e aggiungiamo lo spigolo duale tra il vertice duale corrente e il successivo
            pair<unsigned int, unsigned int> latoDuale_daTrovare = {min(idVerticeDuale_precedente, nuovo_id_verticeDuale), max(idVerticeDuale_precedente, nuovo_id_verticeDuale)};
            auto it_dual_edge = dualEdgeToIdMap.find(latoDuale_daTrovare);

            if (it_dual_edge != dualEdgeToIdMap.end()) {
                latiOrdinati_facciaDuale.push_back(it_dual_edge->second); 
                // lo spigolo viene trovato e quindi aggiunto a latiOrdinati_facciaDuale
            } else {
                cerr << "errore, spigolo duale non trovato per facce (vertici duali) "
                          << idVerticeDuale_precedente << " e " << nuovo_id_verticeDuale << " attorno al vertice originale "
                          << id_verticiOriginali << endl;
                          // Significa che non esiste uno spigolo nella mesh originale che colleghi le due facce currentDualVertexId e
                          // nuovo_id_verticeDuale (che sono i vertici della nostra faccia duale) e che passi attraverso id_verticiOriginali
                break;
            }

            // Se abbiamo raggiunto il punto di partenza (per chiudere il poligono)
            if (nuovo_id_verticeDuale == startid_faccia) {
                if (verticiOrdinati_facciaDuale.size() == facceIncidenti.size()) {
                    // Abbiamo aggiunto tutti i vertici unici, possiamo chiudere il ciclo.
                    break;
                } else {
                    // Siamo tornati all'inizio ma non abbiamo visitato tutte le facce incidenti.
                    cerr << "attenzione, ciclo incompleto per vertice (errore topologico) " << id_verticiOriginali << endl;
                    break;
                }
            }
			
			// Aggiungiamo il prossimo vertice duale (id della faccia originale)
			bool giaAggiunto = false;
            for(unsigned int v : verticiOrdinati_facciaDuale) {
                if (v == nuovo_id_verticeDuale) {
                    giaAggiunto = true;
                    break;
                }
            }
            if (!giaAggiunto) {
                verticiOrdinati_facciaDuale.push_back(nuovo_id_verticeDuale);
            }
            
            // Aggiorniamo idVerticeDuale_precedente per la prossima iterazione
            idVerticeDuale_precedente = nuovo_id_verticeDuale;
			
            // Troviamo il prossimo spigolo originale incidente a id_verticiOriginali (vertice originale attorno al quale stiamo costruendo la faccia duale)
            // che connette currentDualVertexId (la faccia precedente) con nuovo_id_verticeDuale (la faccia successiva).
            // Dobbiamo trovare lo spigolo della faccia `nuovo_id_verticeDuale` che è incidente a `id_verticiOriginali`
            // e che non è `currentid_originale_lato`.
            bool nuovoVerticeDualeCiclo = false;
            const auto& lati_prossimaFaccia = mesh_triangolata.Cell2DsEdges[nuovo_id_verticeDuale]; // spigoli della faccia originale che abbiamoa appena raggiunto
            for (unsigned int lato_prossimaFaccia : lati_prossimaFaccia) {
                if (lato_prossimaFaccia == currentid_originale_lato) continue; // Non prendere lo spigolo da cui siamo venuti

                // Controlliamo se lato_prossimaFaccia contiene id_verticiOriginali
                unsigned int v_e1 = mesh_triangolata.Cell1DsExtrema(lato_prossimaFaccia, 0);
                unsigned int v_e2 = mesh_triangolata.Cell1DsExtrema(lato_prossimaFaccia, 1);

                if (v_e1 == id_verticiOriginali || v_e2 == id_verticiOriginali) {
                    currentid_originale_lato = lato_prossimaFaccia; // Questo è il prossimo spigolo da percorrere
                    currentEdgeKey = {min(v_e1, v_e2), max(v_e1, v_e2)};
                    nuovoVerticeDualeCiclo = true;
                    break;
                }
            }
            if (!nuovoVerticeDualeCiclo) {
                cerr << "Errore, Impossibile trovare il prossimo spigolo per il ciclo attorno al vertice " << id_verticiOriginali << endl;
                break;
            }
        } // Fine del ciclo di costruzione della singola faccia duale

        // Assegniamo la faccia duale costruita
        if (!verticiOrdinati_facciaDuale.empty()) {
            mesh_duale.Cell2DsId[counter_idFacceDuale] = counter_idFacceDuale;
            mesh_duale.Cell2DsVertices[counter_idFacceDuale] = verticiOrdinati_facciaDuale;
            mesh_duale.Cell2DsEdges[counter_idFacceDuale] = latiOrdinati_facciaDuale;
            counter_idFacceDuale++;
        }
	}
}

void proj_sphere(PolyhedralMesh& mesh_triangolata) {
    for (int i = 0; i < mesh_triangolata.Cell0DsCoordinates.cols(); ++i) {
        Vector3d cordinateVertice = mesh_triangolata.Cell0DsCoordinates.col(i);
        double norm = cordinateVertice.norm(); // Equivalente a std::sqrt(cordinateVertice.squaredNorm());
        if (norm < 1e-12) { 
            cerr << "attenzione, Vertice " << i << " troppo vicino all'origine (non proiettatabile)." << endl;
            continue; // Salta la proiezione per questo vertice
        }
        mesh_triangolata.Cell0DsCoordinates.col(i) = cordinateVertice / norm;
    }
}

}