#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <queue>
#include <vector>   
#include <map>                
#include <limits>             
#include <algorithm>          
#include <tuple>              
#include <utility>
#include <Eigen/Dense>          


using namespace std;

using namespace std;
using namespace Eigen;
vector<vector<double>> base_spigoli;
vector<vector<double>> base_vertici;
    int new_id_v = VP.size()+1
    int new_id_s = 1;
    int new_id_f = 1;

vector<double> creatore_vertici(vector<double> A, vector<double> B, vector<double> C, int id){ 
    double diff_x = B[1]-A[1]; 
    double diff_y = B[2]-A[2];
    double diff_z = B[3]-A[3];
    vector<double> nuovo_vertice = { id + 1, C[1] + diff_x, C[2] + diff_y, C[3] + diff_z, 0 };
    id++;
    return nuovo_vertice;
}

void controllo_ordine(vector<vector<double>>& matrice1_faccia, vector<vector<double>>& matrice2_faccia, vector<vector<double>>& matrice3_faccia){   
    if (matrice1_faccia[1][0] == matrice2_faccia[1][0] && matrice1_faccia.back()[0] == matrice3_faccia[1][0]){
        return;
    }
    if (matrice1_faccia[1][0] == matrice2_faccia[1][0]){
        reverse(matrice3_faccia.begin(), matrice3_faccia.end());
        return;
    }
    if (matrice1_faccia.back()[0] == matrice3_faccia[1][0]){
        reverse(matrice2_faccia.begin(), matrice2_faccia.end());
        return;
    }
    reverse(matrice2_faccia.begin(), matrice2_faccia.end());
    reverse(matrice3_faccia.begin(), matrice3_faccia.end());
    return;
}

void mesh1(int new_id_v, int new_id_s, int new_id_f, vector<vector<double>> SP, vector<vector<double>> VP, vector<vector<double>> FP, vector<vector<double>> B, vector<vector<double>> NV, vector<vector<double>> NS, vector<vector<double>> NF){

    NV = VP;
    for (const auto& i : SP){
        vector<int> estremi_spigolo = {i[1],i[2]};
        vector<vector<double>> dati_estremi;
        for (int val : estremi_spigolo) {
            auto riga = find_if(VP.begin(), VP.end(), [val](const vector<int>& riga) {
                return !riga.empty() && riga[0] == val;
            });
            dati_estremi.push_back(*riga);
        }
        base_spigoli.push_back({new_id_s, dati_estremi[0][0], new_id_v});
        base_vertici.push_back({dati_estremi[0][0], dati_estremi[0][1], dati_estremi[0][2], dati_estremi[0][3]});
        base_vertici.push_back({dati_estremi[1][0], dati_estremi[1][1], dati_estremi[1][2], dati_estremi[1][3]});
        new_id_s++;
        for (int j=1, j<b, ++j){
            base_vertici.push_back({new_id_v, min(dati_estremi[0][1], dati_estremi[1][1]) + j*(max(dati_estremi[0][1], dati_estremi[1][1])-min(dati_estremi[0][1], dati_estremi[1][1]))/b,
                          min(dati_estremi[0][2], dati_estremi[1][2]) + j*(max(dati_estremi[0][2], dati_estremi[1][2])-min(dati_estremi[0][2], dati_estremi[1][2]))/b,
                          min(dati_estremi[0][3], dati_estremi[1][3]) + j*(max(dati_estremi[0][3], dati_estremi[1][3])-min(dati_estremi[0][3], dati_estremi[1][3]))/b, i[0]
                        });
            NV.push_back({new_id_v, min(dati_estremi[0][1], dati_estremi[1][1]) + j*(max(dati_estremi[0][1], dati_estremi[1][1])-min(dati_estremi[0][1], dati_estremi[1][1]))/b,
                          min(dati_estremi[0][2], dati_estremi[1][2]) + j*(max(dati_estremi[0][2], dati_estremi[1][2])-min(dati_estremi[0][2], dati_estremi[1][2]))/b,
                          min(dati_estremi[0][3], dati_estremi[1][3]) + j*(max(dati_estremi[0][3], dati_estremi[1][3])-min(dati_estremi[0][3], dati_estremi[1][3]))/b, i[0]
                        });
            if (j<b-1){
                base_spigoli.push_back({new_id_s, new_id_v, new_id_v + 1, i[0]});
                new_id_s++;
            }
            new_id_v++;
        }
        base_spigoli.push_back({new_id_s, new_id_v, dati_estremi[1][0]});
        new_id_s++;
        NS.push_back(base_spigoli);
    }

    for(size_t faccia = 0, i < FP.size(), ++faccia){
        int faccia1 = FP[faccia][4];
        int faccia2 = FP[faccia][5];
        int faccia3 = FP[faccia][6];
        vector<vector<double>> matrice1faccia;
        vector<vector<double>> matrice2faccia;
        vector<vector<double>> matrice3faccia;
        for (const auto& vec : base_vertici) {
            if (!vec.empty() && vec.back() == faccia1) {
                matrice1faccia.push_back(vec); 
            }
            if (!vec.empty() && vec.back() == faccia2) {
                matrice2faccia.push_back(vec); 
            }
            if (!vec.empty() && vec.back() == faccia3) {
                matrice3faccia.push_back(vec); 
            }
        }
        vector<vector<double>> transitore;
        vector<vector<double>> transitore_facce;
        vector<double> transitore_spigolo;
        controllo_ordine(matrice1faccia, matrice2faccia, matrice3faccia);
        for (i = 0, i < b - 1, ++i){
            transitore.push_back(matrice2faccia[i+1]);
            
            for(j = 0, j < b - i - 1, ++j){
                vector<double> nuovo_punto = creatore_vertici(matrice1faccia[j], matrice1faccia[j+1], transitore.back(), new_id_v);
                if (j != b - i - 2){
                    NV.push_back(nuovo_punto);
                }
                else {;
                    nuovo_punto[0] = matrice3faccia[i+1][0];
                }
                if (j>0){
                    transitore_facce.back().push_back(new_id_s);
                    NS.push_back({new_id_s, matrice1faccia[j],transitore.back()[0]});
                    new_id_s++;
                }
                NS.push_back({new_id_s, matrice1faccia[j+1][0], transitore.back()[0]});
                new_id_s;
                NS.push_back({new_id_s, transitore.back()[0], nuovo_punto[0]});
                transitore_facce.push_back({new_id_f, matrice1faccia[j+1][0], transitore.back()[0], nuovo_punto[0], new_id_s});
                new_id_f++;
                transitore_spigolo.push_back(new_id_s);
                new_id_s++;
                NF.push_back({new_id_f, matrice1faccia[j][0], matrice1faccia[j+1][0], transitore.back()[0], new_id_s, new_id_s-1, new_id_s-2});
                new_id_f++;
                transitore.push_back(nuovo_punto);
            }
            NS.push_back({new_id_s, matrice1faccia[matrice1faccia.size() - 2][0], matrice1faccia.back()[0]});
            new_id_s++;
            NS.push_back({new_id_s, matrice1faccia[matrice1faccia.size() - 2][0], transitore.back()[0]});
            transitore_facce.back().push_back(new_id_s);
            new_id_s++;
            NS.push_back({new_id_s, matrice1faccia.back()[0], transitore.back()[0]});
            new_id_s++;
            NF.push_back({new_id_f, matrice1faccia[matrice1faccia.size() - 2][0], matrice1faccia.back()[0], nuovo_punto[0], new_id_s, new_id_s - 1, new_id_s - 2,});
            new_if_f++;

            if (i > 1){    
                for (k = 1, k < b - 1, k++){
                    transitore_facce[k].push_back(transitore_spigolo[k]);
                }
            }
            NF.push_back(transitore_facce);
            matrice1faccia = transitore;
            transitore.clear();
            
            if (i != b - 2){
                transitore_spigolo.clear();
                transitore_facce.clear();
            }
            
        }

        transitore_facce.back().push_back(transitore_spigolo.back());
        NF.push_back(transitore_facce);

        NF.push_back({new_id_f, matrice2faccia.back()[0], matrice2faccia[matrice2faccia.size()-2][0], matrice3faccia[matrice3faccia.size()-2][0]})
        for (const auto& vec : base_spigoli) {
            if (!vec.empty() && vec.back() == faccia2 && vec[2] == matrice2faccia.back()[0] || vec[1] == matrice2faccia.back()[0]){
                NF.back().push_back(vec[0]);
            }
            if (!vec.empty() && vec.back() == faccia3 && vec[2] == matrice3faccia.back()[0] || vec[1] == matrice3faccia.back()[0]){
                NS.back().push_back(vec[0]);
            }
        }    
        NF.back().push_back(transitore_spigolo.back()[0]);
        new_id_f++;       

    }
}



int CheckAddVertices(PolyhedralMesh& poly, const Vector3d& vertex, int& id_vert){
	
	for (unsigned int i = 0; i < poly.Cell0DsId.size(); i++){
		
		if((poly.Cell0DsCoordinates.col(i) - vertex).norm() < 1e-12){
			return i;
		}
	}
	id_vert++;
	poly.Cell0DsId.push_back(id_vert);
	poly.Cell0DsCoordinates.col(id_vert) = vertex; 
	poly.ShortPathCell0Ds[0].push_back(id_vert);
	
	return id_vert;
}


int CheckAddEdges(PolyhedralMesh& poly, const Vector2i& edge, int& id_edge){
	int w0 = edge[0]; 
	int w1 = edge[1]; 
	
	for (unsigned int i = 0; i < poly.Cell1DsId.size(); i++){
		int u0 = poly.Cell1DsExtrema(0,i);
		int u1 = poly.Cell1DsExtrema(1,i);
		
		if((w0 == u0 && w1 == u1)||(w0 == u1 && w1 == u0)){
			return i;//id del lato che verrebbe duplicato
		}
	}
	id_edge++;
	poly.Cell1DsId.push_back(id_edge);
	poly.Cell1DsExtrema.col(id_edge) = edge;
	poly.ShortPathCell1Ds[0].push_back(id_edge);
	
	return id_edge;	
}




void TriangulationTypeII(const PolyhedralMesh& polyOld, PolyhedralMesh& polyNew,const int& n){
	
	int idV_new = -1;// id dei nuovi vertici del poliedro generato dopo la triangolazione
	int idE_new = -1;// id dei nuovi lati del poliedro generato dopo la triangolazione
	int idF_new = 0;// id delle nuove facce del poliedro generato dopo la triangolazione
	
	int numV = polyOld.NumCell0Ds;
	int numE = polyOld.NumCell1Ds;
	int numF = polyOld.NumCell2Ds;
	
	int V = numV + numE*(2*n-1) + numF *(1.5*n*n - 1.5*n +1);
	int E = numE * 2*n + numF*(4.5*n*n + 1.5*n);
	int F = numF * (3*n*n + 3*n);
	
	//predispongo la struttura della PolyhedralMesh andando ad allocare sufficiente memoria per tutte le strutture dati utilizzate per la memorizzazione delle informazioni
	polyNew.NumCell0Ds = V;
	polyNew.Cell0DsId.reserve(V);
	polyNew.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, V);
	
	polyNew.NumCell1Ds = E;
	polyNew.Cell1DsId.reserve(E);
	polyNew.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, E);
	
	polyNew.NumCell2Ds = F;
	polyNew.Cell2DsId.reserve(F);
	polyNew.Cell2DsEdges.reserve(F);
	polyNew.Cell2DsVertices.reserve(F);
	
	polyNew.NumCell3Ds = 1;
	polyNew.Cell3DsId.reserve(1);
	polyNew.Cell3DsEdges.reserve(1);
	polyNew.Cell3DsVertices.reserve(1);
	polyNew.Cell3DsFaces.reserve(1);
		
	for (int f = 0; f < polyOld.NumCell2Ds; f++) {//itero sulle facce del poliedro di base
		
		//mi salvo le coordinate dei vertici della faccia del poliedro di partenza che voglio triangolare
		Vector3d v0 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][0]);
		Vector3d v1 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][1]);
		Vector3d v2 = polyOld.Cell0DsCoordinates.col(polyOld.Cell2DsVertices[f][2]);
		
		//genero i nuovi vertici, assegnando a ognuno una diversa combinazione di coefficienti (coordinate baricentriche)
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		Vector3d new_point(0.0, 0.0, 0.0);
		MatrixXi vertPosition = MatrixXi::Zero(n+1, n+1);
		for(int i=0;i<=n;i++){ //itero per creare tutte le combinazioni di coefficienti possibili 
            for (int j = 0; j<=i;j++){				
                a=double(i-j)/n;
                b=double(j)/n;
                c=1.0-a-b;
					
				new_point = a*v0 + b*v1 + c*v2;
				
				//proietto il nuovo vertice creato sulla superficie di una sfera unitaria centrata nell'origine
				if (new_point.norm() < 1e-16) {
					cerr << "Warning: il vettore considerato ha lunghezza nulla";
					break;}
				new_point.normalize();
				
				//controllo se il punto è già presente
				int id_point = CheckAddVertices(polyNew,new_point,idV_new);
				vertPosition(i,j) = id_point;//inserisce al posto del nuovo id, quello del corrispondente già creato	
			}
		}
		
		//creo una matrice per sfruttare le relazioni con gli indici e collegare i baricentri
		MatrixXi centroidPosition = MatrixXi::Zero(n,2*n-1);
		
		//creo un vettore di vettori che memorizza i 3 vertici dei triangoli "in giù" - utile per relazionarli con i baricentri e creare i triangoli
		vector<vector<int>> verticesToConnect = {};
		verticesToConnect.reserve(n*n-0.5*n*(n+1));
		
		//genero i triangoli (nuovi lati e facce)
		for(int i=0;i < n;i++){
			//indica la colonna del punto in centroidPosition
			int p = 0;
			for (int j = 0; j<=i;j++){
				int vert1 = vertPosition(i,j);
				int vert2 = vertPosition(i+1,j);
				int vert3 = vertPosition(i+1,j+1);
				
				// baricentro del primo triangolo con "punta in su" 
				Vector3d centroid1 = (polyNew.Cell0DsCoordinates.col(vert1) + polyNew.Cell0DsCoordinates.col(vert2) + polyNew.Cell0DsCoordinates.col(vert3)) / 3.0;
				
				if (centroid1.norm() < 1e-16) {
					cerr << "Warning: il vettore considerato ha lunghezza nulla";
					break;}
				centroid1.normalize();
				
				int id_centroid1 = CheckAddVertices(polyNew,centroid1,idV_new);
				centroidPosition(i,p) = id_centroid1;
				p++;
				
				if(j == 0){
					Vector3d midPoint = (polyNew.Cell0DsCoordinates.col(vert1) + polyNew.Cell0DsCoordinates.col(vert2)) / 2.0;
					
					if (midPoint.norm() < 1e-16) {
						cerr << "Warning: il vettore considerato ha lunghezza nulla";
						break;}
					midPoint.normalize();
					
					int id_midPoint = CheckAddVertices(polyNew,midPoint,idV_new);
					
					//triangolo "sopra"
					int id1n = CheckAddEdges(polyNew, {vert1, id_centroid1}, idE_new);
					int id2n = CheckAddEdges(polyNew, {id_centroid1, id_midPoint}, idE_new);
					int id3n = CheckAddEdges(polyNew, {id_midPoint, vert1}, idE_new);
					
					polyNew.Cell2DsVertices.push_back({vert1,id_centroid1,id_midPoint});
					polyNew.Cell2DsEdges.push_back({id1n,id2n,id3n});
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
											
					//triangolo "sotto"
					int id4n = CheckAddEdges(polyNew, {id_midPoint, vert2}, idE_new);
					int id5n = CheckAddEdges(polyNew, {vert2, id_centroid1}, idE_new);
										
					polyNew.Cell2DsVertices.push_back({id_centroid1,id_midPoint,vert2});
					polyNew.Cell2DsEdges.push_back({id2n,id4n,id5n});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
				}
				//creo punto medio su lato destro e creo i due triangoli corrispondenti
				if(i == j){
					Vector3d midPoint = (polyNew.Cell0DsCoordinates.col(vert1) + polyNew.Cell0DsCoordinates.col(vert3)) / 2.0;
										
					if (midPoint.norm() < 1e-16) {
						cerr << "Warning: il vettore considerato ha lunghezza nulla";
						break;}
					midPoint.normalize();
					
					int id_midPoint = CheckAddVertices(polyNew,midPoint,idV_new);

					//triangolo "sopra"
					int id1n = CheckAddEdges(polyNew, {vert1, id_centroid1}, idE_new);
					int id2n = CheckAddEdges(polyNew, {id_centroid1, id_midPoint}, idE_new);
					int id3n = CheckAddEdges(polyNew, {id_midPoint, vert1}, idE_new);
					
					polyNew.Cell2DsVertices.push_back({vert1,id_centroid1,id_midPoint});
					polyNew.Cell2DsEdges.push_back({id1n,id2n,id3n});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
					
					//triangolo "sotto"
					int id4n = CheckAddEdges(polyNew, {id_midPoint, vert3}, idE_new);
					int id5n = CheckAddEdges(polyNew, {vert3, id_centroid1}, idE_new);
					
					
					polyNew.Cell2DsVertices.push_back({id_centroid1,id_midPoint,vert3});
					polyNew.Cell2DsEdges.push_back({id2n,id4n,id5n});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
				}
				//creo punto medio su lato inferiore e creo i due triangoli corrispondenti
				if(i == n-1){
					Vector3d midPoint = (polyNew.Cell0DsCoordinates.col(vert2) + polyNew.Cell0DsCoordinates.col(vert3)) / 2.0;
					
					if (midPoint.norm() < 1e-16) {
						cerr << "Warning: il vettore considerato ha lunghezza nulla";
						break;}
					midPoint.normalize();
					
					int id_midPoint = CheckAddVertices(polyNew,midPoint,idV_new);
					
					//triangolo "destra"
					int id1n = CheckAddEdges(polyNew, {vert2, id_centroid1}, idE_new);
					int id2n = CheckAddEdges(polyNew, {id_centroid1, id_midPoint}, idE_new);
					int id3n = CheckAddEdges(polyNew, {id_midPoint, vert2}, idE_new);
					
					polyNew.Cell2DsVertices.push_back({vert2,id_centroid1,id_midPoint});
					polyNew.Cell2DsEdges.push_back({id1n,id2n,id3n});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
					
					//triangolo "sinistra"
					int id4n = CheckAddEdges(polyNew, {id_midPoint, vert3}, idE_new);
					int id5n = CheckAddEdges(polyNew, {vert3, id_centroid1}, idE_new);
					
					polyNew.Cell2DsVertices.push_back({id_centroid1,id_midPoint,vert3});
					polyNew.Cell2DsEdges.push_back({id2n,id4n,id5n});
					
					polyNew.Cell2DsId.push_back(idF_new);
					idF_new++;
				}
				if(i!=j){
					int vert4 = vertPosition(i,j+1);
					
					verticesToConnect.push_back({vert1,vert3,vert4});
					
					// baricentro del primo triangolo con "punta in giù" 
					Vector3d centroid2 = (polyNew.Cell0DsCoordinates.col(vert1) + polyNew.Cell0DsCoordinates.col(vert3) + polyNew.Cell0DsCoordinates.col(vert4)) / 3.0;
					
					if (centroid2.norm() < 1e-16) {
						cerr << "Warning: il vettore considerato ha lunghezza nulla";
						break;}
					centroid2.normalize();
					
					int id_centroid2 = CheckAddVertices(polyNew,centroid2,idV_new);
					
					centroidPosition(i,p) = id_centroid2;
					p++;
				}					
			}
		}
		
		int k = 0;
		for(int i = 1;i<n;i++){
			for(int j=1;j<=2*i;j+=2){
				int id_c = centroidPosition(i,j); //prendo l'id del baricentro che voglio collegare
				
				//memorizzo vertici del triangolo considerato
				int vx0 = verticesToConnect[k][0];
				int vx1 = verticesToConnect[k][1];
				int vx2 = verticesToConnect[k][2];
				
				//memorizzo i baricentri circostanti da collegare
				int c_up = centroidPosition(i-1,j-1);
				int c_left = centroidPosition(i,j-1);
				int c_right = centroidPosition(i,j+1);
				
				//creo i nuovi lati dei triangoli intorno al baricentro considerato
				int i1 = CheckAddEdges(polyNew, {vx0, c_up}, idE_new);
				int i2 = CheckAddEdges(polyNew, {c_up, id_c}, idE_new);
				int i3 = CheckAddEdges(polyNew, {id_c, vx0}, idE_new);
				int i4 = CheckAddEdges(polyNew, {id_c, vx2}, idE_new);
				int i5 = CheckAddEdges(polyNew, {vx2, c_up}, idE_new);
				int i6 = CheckAddEdges(polyNew, {vx2, c_right}, idE_new);
				int i7 = CheckAddEdges(polyNew, {c_right, id_c}, idE_new);
				int i8 = CheckAddEdges(polyNew, {id_c, vx1}, idE_new);
				int i9 = CheckAddEdges(polyNew, {vx1, c_right}, idE_new);
				int i10 = CheckAddEdges(polyNew, {vx1, c_left}, idE_new);				
				int i11 = CheckAddEdges(polyNew, {c_left, id_c}, idE_new);
				int i12 = CheckAddEdges(polyNew, {vx0, c_left}, idE_new);
				
				//creo i 6 nuovi triangoli intorno al baricentro
				//triangolo 1
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({vx0,c_up,id_c});
				polyNew.Cell2DsEdges.push_back({i1,i2,i3});
				idF_new++;
				
				//triangolo 2
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({c_up,id_c,vx2});
				polyNew.Cell2DsEdges.push_back({i2,i4,i5});
				idF_new++;
				
				//triangolo 3
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({id_c,vx2,c_right});
				polyNew.Cell2DsEdges.push_back({i4,i6,i7});
				idF_new++;
				
				//triangolo 4
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({c_right,id_c,vx1});
				polyNew.Cell2DsEdges.push_back({i7,i8,i9});
				idF_new++;
				
				//triangolo 5
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({id_c,vx1,c_left});
				polyNew.Cell2DsEdges.push_back({i8,i10,i11});
				idF_new++;
				
				//triangolo 6
				polyNew.Cell2DsId.push_back(idF_new);
				polyNew.Cell2DsVertices.push_back({id_c,vx0,c_left});
				polyNew.Cell2DsEdges.push_back({i3,i12,i11});
				idF_new++;
				k++;
			}
		}
	}
	//inserisco nella PolyhedralMesh le informazioni su Cell3Ds
	polyNew.Cell3DsId.push_back(0); 
	polyNew.Cell3DsVertices.push_back(polyNew.Cell0DsId);
	polyNew.Cell3DsEdges.push_back(polyNew.Cell1DsId);
	polyNew.Cell3DsFaces.push_back(polyNew.Cell2DsId);
}

// IMPLEMENTAZIONE DEL COSTRUTTORE DI SHORTESTPATHRESULT
ShortestPathResult::ShortestPathResult(unsigned int nEdges, double len)
    : numEdges(nEdges), totalLength(len)
{}

// Calcolo la distanza euclidea tra due punti.
double calculateDistanceById(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2) {
    
    VectorXd A = mesh.M0D[0].size(id1);
    VectorXd B = mesh.M0D[0].size(id2);
	double Norma = (A - B).norm(); // Calcola la norma (distanza euclidea)

    return Norma; 
}


// Funzione per calcolare la matrice di adiacenza
MatrixXi CalcoloMatriceAdiacenza(const PolyhedralMesh& mesh) {
	
    const unsigned int numVertices = mesh.M0D[0].size();
	MatrixXi Matr_Adiacenza = MatrixXi::Zero(numVertices, numVertices);

    // Itera su tutti i lati (Cell1Ds) della mesh
    for (unsigned int i = 0; i < mesh.M1D.size(); ++i) {

        unsigned int v1 = mesh.M1D(i, 1);
        unsigned int v2 = mesh.M1D(i, 2);
		
        Matr_Adiacenza(v1, v2) = 1;
        Matr_Adiacenza(v2, v1) = 1;
    }

    return Matr_Adiacenza;
}

ShortestPathResult findShortestPathDijkstra(
    PolyhedralMesh& mesh,
    const MatrixXi& Matr_Adiacenza,
    const unsigned int nodo_i,
    const unsigned int nodo_f ) {
		
	const unsigned int numVertici = mesh.M0D[0].size(); // Numero totale di vertici
    const unsigned int numLatiMesh = mesh.M1D[0].size(); // Numero totale di lati nella mesh
    
    if (nodo_i >= numVertici) {
        cerr << "Errore: il vertice di partenza (" << nodo_i << ") è fuori dal range [0, " << numVertici - 1 << "]." << endl;
        // Restituisce un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }
    if (nodo_f >= numVertici) {
        cerr << "Errore: il vertice di arrivo (" << nodo_f << ") è fuori dal range [0, " << numVertici - 1 << "]." << endl;
        // Restituisce un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }

	// Inizializza il risultato
    ShortestPathResult result(0, 0.0);

    // Caso in cui partenza e arrivo sono uguali
    if (nodo_i == nodo_f) {
        cout << "Partenza e arrivo coincidono quindi il cammino è nullo." << endl;
        mesh.M0DMarker.assign(numVertici, 0); 
        mesh.M1DMarker.assign(numLatiMesh, 0);
        mesh.M0DMarker[nodo_i] = 1; // Assegna il marker al vertice sulla mesh
        return result;
    }

    // Mappa per collegare una coppia di vertici (tramite INDICI)
    // all'ID del lato e alla sua lunghezza.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> MappaLati;
    
	for (unsigned int i = 0; i < numLatiMesh; ++i) {
        unsigned int v1 = mesh.M1D(i, 1);
        unsigned int v2 = mesh.M1D(i, 2);

        double lunghezza = calculateDistanceById(mesh, v1, v2);

        // Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
        MappaLati[{min(v1, v2), max(v1, v2)}] = {mesh.M1D[i][0], lunghezza};
    }

    // Variabili Dijkstra
	
    // dist[i] = distanza minima conosciuta dal vertice di partenza a i
    vector<double> dist(numVertici, numeric_limits<double>::infinity());
	
    // predVertex[i] = indice del vertice precedente nel cammino minimo a i
    vector<unsigned int> predVertex(numVertici, -1); // Usiamo -1 = nessun predecessore
	
    // predEdge[i] = ID del Cell1D usato per raggiungere i dal suo predecessore
    vector<unsigned int> predEdge(numVertici, -1);

    // visited[i] = true se il cammino più breve a 'i' è stato finalizzato
    vector<bool> visitato(numVertici, false); 

    // Coda di priorità: memorizza coppie {distanza, indice_vertice}
    // std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
    using QueueElementi = pair<double, unsigned int>;
    priority_queue<QueueElementi, vector<QueueElementi>, greater<>> pq;
	
	// priority_queue<QueueElem, vector<QueueElem>, greater<QueueElem>> pq;

    // Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
    dist[nodo_i] = 0.0;
    pq.push({0.0, nodo_i});

    // algoritmo Dijkstra
    while (!pq.empty()) {
        // Estrai il vertice 'u' con la distanza minima corrente dalla coda
		auto [current_distance, u] = pq.top(); // accedo all'elemento
        pq.pop(); // rimuovo l'elemento

        // Se il vertice è già stato visitato, abbiamo già trovato
        // il cammino più breve per esso, quindi ignoriamo.
        if (visitato[u]) {
            continue;
        }
        visitato[u] = true; // Marca il vertice come visitato/finalizzato

        // Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
        if (u == nodo_f) {
            break;
        }

        // Itera su tutti i possibili vicini 'v' di 'u' usando la matrice di adiacenza
        for (unsigned int v = 0; v < numVertici; ++v) {
            // Se c'è un lato tra u e v (adjMatrix(u, v) == 1)
            if (Matr_Adiacenza(u, v) == 1) {
                // Recupera le informazioni sul lato (ID reale e lunghezza)
                auto it_edge = MappaLati.find({min(u, v), max(u, v)});
                if (it_edge == MappaLati.end()) {
                    cerr << "Attenzione: Lato tra indici (" << u << "," << v << ") presente in matrice di adiacenza ma non trovato in MappaLati.\n";
                    continue;
                }
                double peso = it_edge->second.second; // Peso del lato = lunghezza euclidea

                // Operazione di "rilassamento":
                // Se la distanza calcolata a 'v' passando per 'u' è minore della distanza attuale di 'v'
                if (dist[u] + peso < dist[v]) {
                    dist[v] = dist[u] + peso; // Aggiorna la distanza minima per 'v'
                    predVertex[v] = u;          // Imposta 'u' come predecessore di 'v'
                    predEdge[v] = it_edge->second.first; // Memorizza l'ID reale del lato usato
                    pq.push({dist[v], v});      // Inserisci 'v' nella coda di priorità con la nuova distanza
                }
            }
        }
    }
	
    // Se la distanza al vertice di arrivo è ancora infinito, significa che non è stato trovato alcun cammino
    if (dist[nodo_f] == numeric_limits<double>::infinity()) {
        cout << "Non è possibile trovare un cammino tra il vertice " << nodo_i
                  << " e il vertice " << nodo_f << endl;
        mesh.M0DMarker.assign(numVertici, 0);
        mesh.M1DMarker.assign(numLatiMesh, 0);
		return result;
    }
	
    // Ricostruiamo il cammino dal vertice di arrivo al vertice di partenza
	
	// Inizializzo i marker della mesh a 0
    mesh.M0DMarker.assign(numVertici, 0); // Inizializza i marker a 0
    mesh.M1DMarker.assign(numLatiMesh, 0);
	
	result.DistanzaTot = dist[nodo_f]; // distanza totale

    unsigned int id_corrente = nodo_f; // Partiamo dall'indice del vertice di arrivo
	
	while (id_corrente != nodo_i) {
        // Marca il vertice corrente come parte del cammino minimo
        mesh.M0DMarker[id_corrente] = 1;

        unsigned int prec_vert_id = predVertex[id_corrente]; // Indice del vertice precedente nel cammino
        unsigned int lato_used_id = predEdge[id_corrente];  // ID del lato usato per raggiungere id_corrente

        mesh.M1DMarker[lato_used_id] = 1;
        result.numEdges++;
        
        current_idx = prec_vert_id; // Spostati al vertice precedente e continua la ricostruzione
    }
	
    // Marca anche il vertice di partenza come parte del cammino
    mesh.M0DMarker[nodo_i] = 1;	
	
	return result;
}