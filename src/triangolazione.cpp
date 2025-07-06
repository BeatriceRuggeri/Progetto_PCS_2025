#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void edge_correction(
		const int a, const int b,
		PolyhedralMesh& meshTriangolazione,
		unsigned int& k2,
		const unsigned int k3)
	{
		bool found = false;
	
		for (unsigned int i = 0; i <= k2; ++i) {  // Scorri solo fino all'ultimo edge inserito
			if ((meshTriangolazione.Cell1DsExtrema(i, 0) == a && meshTriangolazione.Cell1DsExtrema(i, 1) == b) ||
				(meshTriangolazione.Cell1DsExtrema(i, 0) == b && meshTriangolazione.Cell1DsExtrema(i, 1) == a)) {
				
				// Edge già presente ⇒ usalo
				meshTriangolazione.Cell2DsEdges[k3].push_back(i);
				found = true;
				break;
			}
		}
		
		if (!found) {
			// Edge non esiste ⇒ lo creiamo
			meshTriangolazione.Cell1DsExtrema.row(k2) << a, b;
			meshTriangolazione.Cell1DsId[k2] = k2;
		
			const auto& c = meshTriangolazione.Cell0DsFlag[a];
			const auto& d = meshTriangolazione.Cell0DsFlag[b];
		
			bool in_comune = false;
			for (size_t s = 0; s < c.size(); ++s) {
				for (size_t t = 0; t < d.size(); ++t) {
					if (c[s] == d[t]) {
						meshTriangolazione.Cell1DsFlag[k2] = c[s];
						in_comune = true;
						break;
					}
				}
				if (in_comune) break;
			}
			
			if (!in_comune) {
				meshTriangolazione.Cell1DsFlag[k2] = numeric_limits<unsigned int>::max();
			}
		
			meshTriangolazione.Cell2DsEdges[k3].push_back(k2);
			++k2;
		}
	}
	

	void gestione_input(PolyhedralMesh& meshT, const vector<int>& dim)
	{
		meshT.Cell0DsId.resize(dim[0]);
        meshT.Cell0DsCoordinates = MatrixXd::Zero(3, dim[0]);
        meshT.Cell0DsFlag.resize(dim[0]);

        meshT.Cell1DsId.resize(dim[1]);
        meshT.Cell1DsExtrema = MatrixXi::Zero(dim[1], 2);
        meshT.Cell1DsFlag.resize(dim[1]);

        meshT.Cell2DsId.resize(dim[2]);
        meshT.Cell2DsVertices.resize(dim[2]);
        meshT.Cell2DsEdges.resize(dim[2]);
		return;

	}

	// Funzione principale per la triangolazione e salvataggio nella mesh
	void tri_build_I(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolazione, const unsigned int b, const unsigned int c, const vector<int>& dimensioneDuplicati) {
		
		unsigned int liv_di_suddivisione = b+c;
		
		gestione_input(meshTriangolazione, dimensioneDuplicati);

		unsigned int k1=0;
		unsigned int k2=0;
		unsigned int k3=0;
		
		for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& face = mesh.Cell2DsVertices[faceId];
			Vector3d V0 = mesh.Cell0DsCoordinates.col(face[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(face[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(face[2]);
	
			vector<vector<int>> vertexGrid;
	
			for (unsigned int i = 0; i <= liv_di_suddivisione; ++i) {
				vector<int> row;
				Vector3d start = ((double)i / liv_di_suddivisione) * V1 + ((double)(liv_di_suddivisione - i) / liv_di_suddivisione) * V0;
				Vector3d end   = ((double)i / liv_di_suddivisione) * V2 + ((double)(liv_di_suddivisione - i) / liv_di_suddivisione) * V0;
	
				for (unsigned int j = 0; j <= i; ++j) {
					Vector3d point;
					if (i == 0) {
						point = V0;
					} else {
						point = ((double)j / i) * end + ((double)(i - j) / i) * start;
					}
					
					meshTriangolazione.Cell0DsCoordinates.col(k1) = point;
					meshTriangolazione.Cell0DsId[k1]=k1;

					if (i == 0) {
					meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][2]};
					}
					else if (i == liv_di_suddivisione) {
						if (j == 0)
							meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][1]};
						else if (j == liv_di_suddivisione)
							meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1], mesh.Cell2DsEdges[faceId][2]};
						else
							meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][1]};
					}
					else if (j == 0) {
						meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][0]};
					}
					else if (j == i) {
						meshTriangolazione.Cell0DsFlag[k1] = {mesh.Cell2DsEdges[faceId][2]};
					}
					else {
						meshTriangolazione.Cell0DsFlag[k1] = {numeric_limits<unsigned int>::max()};
					}		
	
					row.push_back(k1);
					k1++;
				}
				vertexGrid.push_back(row);
			}
			for (unsigned int i = 0; i < liv_di_suddivisione; ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					unsigned int v1 = vertexGrid[i][j];
					unsigned int v2 = vertexGrid[i + 1][j];
					unsigned int v3 = vertexGrid[i + 1][j + 1];
		
					vector<unsigned int> verts1 = {v1, v2, v3};
					meshTriangolazione.Cell2DsVertices[k3]=verts1;
						
					meshTriangolazione.Cell2DsId[k3]=k3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts1[e];
						int b = verts1[(e + 1) % 3];
						edge_correction(a, b, meshTriangolazione, k2, k3);
					}
	
					k3++;
							
					unsigned int v4 = vertexGrid[i][j];
					unsigned int v5 = vertexGrid[i + 1][j + 1];
					unsigned int v6 = vertexGrid[i][j + 1];
		
					vector<unsigned int> verts2 = {v4, v5, v6};

					meshTriangolazione.Cell2DsVertices[k3]=verts2;
					meshTriangolazione.Cell2DsId[k3]=k3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts2[e];
						int b = verts2[(e + 1) % 3];
						edge_correction(a, b, meshTriangolazione, k2, k3);
					}
					k3++;
				}
		
				// Triangolo finale in basso a sinistra
				unsigned int v1 = vertexGrid[i][i];
				unsigned int v2 = vertexGrid[i + 1][i];
				unsigned int v3 = vertexGrid[i + 1][i + 1];
	
				vector<unsigned int> verts = {v1, v2, v3};
				meshTriangolazione.Cell2DsVertices[k3]=verts;

				meshTriangolazione.Cell2DsId[k3]=k3;
					
				for (unsigned int e = 0; e < 3; ++e) {
					int a = verts[e];
					int b = verts[(e + 1) % 3];
					edge_correction(a, b, meshTriangolazione, k2, k3);
				}
				k3++;
		}
			
		}

	}
	
	// Funzione per trovare la forma normalizzata ciclicamente
	// cioè, ruotare la sequenza in modo che inizi con l'elemento più piccolo


	vector<unsigned int> norm_cyclic(const vector<unsigned int>& current_edges) {
		if (current_edges.empty() || current_edges.size() == 1) {
			return current_edges; // Non c'è nulla da normalizzare per 0 o 1 elemento
		}
		vector<unsigned int> temp_edges = current_edges; // Lavora su una copia
		auto min_it = min_element(temp_edges.begin(), temp_edges.end());
		rotate(temp_edges.begin(), min_it, temp_edges.end());
		return temp_edges;
	}

	vector<unsigned int> normalization_L(const vector<unsigned int>& face_edges) {
		if (face_edges.empty() || face_edges.size() == 1) {
			return face_edges; // Non c'è nulla da normalizzare per 0 o 1 elemento
		}
	
		// Ottieni la forma normalizzata ciclicamente della sequenza originale
		vector<unsigned int> normalizzazione_originale = norm_cyclic(face_edges);
	
		// Ottieni la sequenza inversa dell'originale
		vector<unsigned int> reversed_edges = face_edges;
		reverse(reversed_edges.begin(), reversed_edges.end()); // Inverte la copia
	
		// Ottieni la forma normalizzata ciclicamente anche della sequenza inversa
		vector<unsigned int> normalizzazione_contrario = norm_cyclic(reversed_edges);
	
		// Confronta le due forme normalizzate e restituisci la minore
		// confronta gli elementi uno a uno
		if (normalizzazione_originale < normalizzazione_contrario) {
			return normalizzazione_originale;
		} else {
			return normalizzazione_contrario;
		}
	}

	void face_F(const vector<unsigned int>& vertici_nuova_faccia,
							   const vector<unsigned int>& spigoli_nuova_faccia,
							   PolyhedralMesh& meshTriangolata,
							   unsigned int& k3)
	{
		// Normalizza la sequenza degli spigoli
		vector<unsigned int> nuovi_spigoli_normalizzati = normalization_L(spigoli_nuova_faccia);
	
		// Itero attraverso le facce esistenti per trovare un duplicato
		bool found = false;
		// mi serve perchè se no incrementerei k3 tutte le volte anche quando la faccia esiste già
		for (unsigned int i = 0; i < meshTriangolata.Cell2DsId.size(); ++i) {
			// Normalizza gli spigoli della faccia corrente per il confronto
			vector<unsigned int> spigoli_originali_normalizzati = normalization_L(meshTriangolata.Cell2DsEdges[i]);
	
			// Confronta le sequenze normalizzate. Se sono uguali, la faccia esiste già.
			if (nuovi_spigoli_normalizzati == spigoli_originali_normalizzati) {
				found = true;
				break;
			}
		}
	
		// Se il ciclo è terminato e non siamo usciti dalla funzione,
		// significa che la faccia non è stata trovata. La aggiungo.
		if (!found){
			meshTriangolata.Cell2DsVertices[k3] = vertici_nuova_faccia;
			meshTriangolata.Cell2DsEdges[k3] = spigoli_nuova_faccia; // Salva la versione non normalizzata
			meshTriangolata.Cell2DsId[k3] = k3; // Assegna il nuovo ID
		
			k3++; // Restituisci il vecchio valore di k3, poi incrementalo per la prossima faccia
		}
	}

	
	unsigned int vertex_correction(const Vector3d& coord, PolyhedralMesh& meshTriangolata, unsigned int& k1)
	{
		double tol = 1e-12;
		for (unsigned int i = 0; i < meshTriangolata.Cell0DsCoordinates.cols(); i++) {
			if ((meshTriangolata.Cell0DsCoordinates.col(i) - coord).norm() < tol) {
				return i;
			}
		}

		// se il ciclo è terminato e non siamo usciti dalla funzione vuole dire che il vertice non è stato trovato
		meshTriangolata.Cell0DsCoordinates.col(k1) = coord;
		meshTriangolata.Cell0DsId[k1] = k1;
		return k1++;
	}
	
	
	unsigned int edge_correction_II(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangolata, unsigned int& k2)
	{
		// Controllo se l'edge (a,b) o (b,a) esiste già
		for (unsigned int i = 0; i < meshTriangolata.Cell1DsExtrema.rows(); ++i) {
			if ((meshTriangolata.Cell1DsExtrema(i, 0) == static_cast<int>(a) && meshTriangolata.Cell1DsExtrema(i, 1) == static_cast<int>(b)) ||
    			(meshTriangolata.Cell1DsExtrema(i, 0) == static_cast<int>(b) && meshTriangolata.Cell1DsExtrema(i, 1) == static_cast<int>(a))) {

				// Edge già presente ⇒ associarlo alla faccia
				return i; // restituisci l'ID dell'edge
			}
		}

		// Edge non presente ⇒ lo aggiungo
		meshTriangolata.Cell1DsExtrema(k2, 0) = a;
		meshTriangolata.Cell1DsExtrema(k2, 1) = b;
		meshTriangolata.Cell1DsId[k2] = k2;
		return k2++; // restituisci il vecchio valore di k2, poi incrementalo
	}
	
	
	// trova la faccia adiacente a face tramite edge, e poi restituire il baricentro di questa faccia adiacente
	Vector3d baricenter_F(const PolyhedralMesh& meshTriangolata, const unsigned int edgeId, const unsigned int idFacciaCorrente, map<pair<unsigned int, unsigned int>, vector<unsigned int>> mappa_verticeAfaccia) {

		// Ottieni i vertici che compongono il lato
		unsigned int v1_id = meshTriangolata.Cell1DsExtrema(edgeId, 0);
		unsigned int v2_id = meshTriangolata.Cell1DsExtrema(edgeId, 1);
		pair<unsigned int, unsigned int> edgeKey = {min(v1_id, v2_id), max(v1_id, v2_id)};

		// Trova le facce associate a questo lato
		auto it = mappa_verticeAfaccia.find(edgeKey);
		if (it == mappa_verticeAfaccia.end()) {
			cerr << "Errore: Edge " << edgeId << " non trovato nella mappa delle facce adiacenti." << endl;
			return Vector3d::Zero(); // Nessuna faccia associata a questo edge
		}

		const vector<unsigned int>& facce_confinanti = it->second;

		if (facce_confinanti.size() == 2) {
			// Ci sono due facce. Una è idFacciaCorrente, l'altra è quella che cerchiamo.
			unsigned int faceId1 = facce_confinanti[0];
			unsigned int faceId2 = facce_confinanti[1];

			unsigned int targetFaceId;
			if (faceId1 == idFacciaCorrente) {
				targetFaceId = faceId2;
			} else if (faceId2 == idFacciaCorrente) {
				targetFaceId = faceId1;
			} else {
				// Questo caso significa che idFacciaCorrente non è una delle facce che condividono l'edge dato.
				cerr << "Errore: Faccia corrente " << idFacciaCorrente << " non condivide l'edge " << edgeId << "." << endl;
				return Vector3d::Zero();
			}

			// Calcola e restituisci il baricentro della faccia adiacente.
			return baricenter_Fetch(meshTriangolata, targetFaceId);

		} else if (facce_confinanti.size() == 1) {
			cerr << "Warning: Edge " << edgeId << " è un bordo della mesh. Nessuna faccia adiacente trovata." << endl;
			return Vector3d::Zero(); 
		} else {
			cerr << "Errore: Edge " << edgeId << " è condiviso da " << facce_confinanti.size() << " facce" << endl;
			return Vector3d::Zero();
		}
    }

    void tri_build_II(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangolata, const vector<int>& dimension, map<pair<unsigned int, unsigned int>, vector<unsigned int>> mappa_verticeAfaccia) {

        gestione_input(meshTriangolata, dimension);

        unsigned int k1 = 0; // Contatore nodi globali (0D)
        unsigned int k2 = 0; // Contatore edge globali (1D)
        unsigned int k3 = 0; // Contatore facce/triangoli globali (2D)
        
        for (unsigned faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& faceVertices = mesh.Cell2DsVertices[faceId]; // id dei vertici della faccia originale
				
			Vector3d V0 = mesh.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(faceVertices[2]);
			
			Vector3d baricentro = ((V0 + V1 + V2)/3.0);
			
			const auto& faceEdges = mesh.Cell2DsEdges[faceId]; // id dei lati della faccia originale

			for (unsigned int e = 0; e <3; e++){
				
				if (!mesh.Cell1DsOriginalFlag[faceEdges[e]]) { // lato di bordo
					Vector3d punto_medio = (mesh.Cell0DsCoordinates.col(faceVertices[e])+mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]))/2.0;

					// TRIANGOLO SOPRA IL PUNTO MEDIO

					unsigned int vertex = vertex_correction(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangolata, k1);
					unsigned int medium = vertex_correction(punto_medio, meshTriangolata, k1);
					unsigned int bar = vertex_correction(baricentro, meshTriangolata, k1);
					
					vector<unsigned int> vertici_nuova_faccia = {vertex, medium, bar};
					
					// lato vertice - punto medio
					unsigned int edgeA = edge_correction_II(vertex, medium, meshTriangolata, k2);
					
					// lato punto medio- baricentro
					unsigned int edgeB = edge_correction_II(medium, bar, meshTriangolata, k2);
					
					// lato baricentro - vertice
					unsigned int edgeC = edge_correction_II(bar, vertex, meshTriangolata, k2);
					
					vector<unsigned int> spigoli_nuova_faccia = {edgeA, edgeB, edgeC};

					face_F(vertici_nuova_faccia, spigoli_nuova_faccia, meshTriangolata, k3);


					// TRIANGOLO SOTTO IL PUNTO MEDIO
					
					unsigned int vertex1 = vertex_correction(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangolata, k1);
					vector<unsigned int> vertici_nuova_faccia_1 = {vertex1, medium, bar};
					
					// lato vertice - punto medio
					unsigned int edgeD = edge_correction_II(vertex1, medium, meshTriangolata, k2);
					
					// lato baricentro - vertice
					unsigned int edgeE = edge_correction_II(bar, vertex1, meshTriangolata, k2);
					
					vector<unsigned int> spigoli_nuova_faccia_1 = {edgeD, edgeB, edgeE};
					
					face_F(vertici_nuova_faccia_1, spigoli_nuova_faccia_1, meshTriangolata, k3);
					
				
				} else {
					
				// TRIANGOLO A SINISTRA
				Vector3d bar2_coord = baricenter_F(mesh, faceEdges[e], faceId, mappa_verticeAfaccia);
				unsigned int vertex2 = vertex_correction(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangolata, k1);
				unsigned int bar1 = vertex_correction(baricentro, meshTriangolata, k1);
				unsigned int bar2 = vertex_correction(bar2_coord, meshTriangolata, k1);
				
				vector<unsigned int> vertici_nuova_faccia2 = {vertex2, bar1, bar2};
				
				// lato vertice - bar1
				unsigned int edgeF = edge_correction_II(vertex2, bar1, meshTriangolata, k2);
					
				// lato bar1- bar2
				unsigned int edgeG = edge_correction_II(bar1, bar2, meshTriangolata, k2);
					
				// lato bar2 - vertice
				unsigned int edgeH = edge_correction_II(bar2, vertex2, meshTriangolata, k2);
					
				vector<unsigned int> new_face_edges2 = {edgeF, edgeG, edgeH};
				
				face_F(vertici_nuova_faccia2, new_face_edges2, meshTriangolata, k3);
					
				
				// TRIANGOLO A DESTRA
				unsigned int vertex3 = vertex_correction(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangolata, k1);
				vector<unsigned int> vertici_nuova_faccia3 = {vertex3, bar1, bar2};
				
				// lato vertice - bar1
				unsigned int edgeI = edge_correction_II(vertex3, bar1, meshTriangolata, k2);
					
				// lato bar2 - vertice
				unsigned int edgeL = edge_correction_II(bar2, vertex3, meshTriangolata, k2);
					
				vector<unsigned int> new_face_edges3 = {edgeI, edgeG, edgeL};
				
				face_F(vertici_nuova_faccia3, new_face_edges3, meshTriangolata, k3);	
				
				}
			}
		}
	}
	

}