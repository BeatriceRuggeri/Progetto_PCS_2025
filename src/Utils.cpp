#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;


void gen_tetraedro(PolyhedralMesh& mesh){
    mesh.M0D = {
			{0, 1.0 / sqrt(3.0), 1.0 / sqrt(3.0), 1.0 / sqrt(3.0)}, 
			{1, 1.0 / sqrt(3.0), -1.0 / sqrt(3.0), -1.0 / sqrt(3.0)}, 
			{2, -1.0 / sqrt(3.0), 1.0 / sqrt(3.0), -1.0 / sqrt(3.0)}, 
			{3, -1.0 / sqrt(3.0), -1.0 / sqrt(3.0), 1.0 / sqrt(3.0)} 
		};
    mesh.M1D = {
			{0, 0, 1}, 
			{1, 0, 2}, 
			{2, 0, 3}, 
			{3, 1, 2}, 
			{4, 1, 3}, 
			{5, 2, 3} 
		};
    mesh.M2D = {
			{0, 0, 2, 1, 1, 3, 0}, 
			{1, 0, 1, 3, 0, 4, 2}, 
			{2, 0, 3, 2, 2, 1, 5}, 
			{3, 1, 2, 3, 5, 3, 4} 
		};
    return;
}

void gen_ottaedro(PolyhedralMesh& mesh){
    mesh.M0D = {
			{0, 1.0, 0.0, 0.0}, 
			{1, -1.0, 0.0, 0.0}, 
			{2, 0.0, 1.0, 0.0},
			{3, 0.0, -1.0, 0.0},
			{4, 0.0, 0.0, 1.0}, 
			{5, 0.0, 0.0, -1.0} 
		};
    mesh.M1D = {
			{0, 0, 2}, 
			{1, 0, 3}, 
			{2, 0, 4}, 
			{3, 0, 5}, 
			{4, 1, 2}, 
			{5, 1, 3},
			{6, 1, 4}, 
			{7, 1, 5}, 
			{8, 2, 4}, 
			{9, 2, 5}, 
			{10, 3, 4}, 
			{11, 3, 5} 
		};
    mesh.M2D = {
			{0, 0, 2, 4, 0, 8, 2}, 
			{1, 0, 3, 4, 2, 1, 10},  
			{2, 1, 3, 4, 10, 5, 6}, 
			{3, 1, 2, 4, 6, 8, 4}, 
			{4, 1, 2, 5, 4, 7, 9}, 
			{5, 2, 0, 5, 9, 0, 3}, 
			{6, 0, 3, 5, 3, 1, 11}, 
			{7, 3, 1, 5, 11, 5, 7} 
		};
	return;	
}

void gen_icosaedro(PolyhedralMesh& mesh){
    mesh.M0D =  {
			{0, 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, 
			{1, 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, 
			{2, 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))}, 
			{3, 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))}, 
			{4, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0}, 
			{5, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0}, 
			{6, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0}, 
			{7, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0}, 
			{8, (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, 
			{9, (1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1), 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))},
			{10, -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0, 1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)}, 
			{11, -((1.0 + sqrt(5.0)) / 2.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1)), 0.0, -(1.0/sqrt((1.0 + sqrt(5.0)) / 2.0^2 + 1))} 
		};
    mesh.M1D = {
			{0, 0, 1},  
			{1, 0, 4},  
			{2, 0, 6}, 
			{3, 0, 8},  
			{4, 0, 10},  
			{5, 1, 5}, 
			{6, 1, 8},  
			{7, 1, 7}, 
			{8, 1, 10}, 
			{9, 2, 3}, 
			{10, 2, 4}, 
			{11, 2, 6}, 
			{12, 2, 9}, 
			{13, 2, 11}, 
			{14, 3, 5}, 
			{15, 3, 7}, 
			{16, 3, 9}, 
			{17, 3, 11},
			{18, 4, 6}, 
			{19, 4, 8}, 
			{20, 4, 9}, 
			{21, 5, 7}, 
			{22, 5, 8}, 
			{23, 5, 9}, 
			{24, 6, 10}, 
			{25, 6, 11}, 
			{26, 7, 10}, 
			{27, 7, 11}, 
			{28, 8, 9}, 
			{29, 10, 11} 
		};
    mesh.M2D = {
			{0, 0, 1, 8, 3, 6, 0},               
			{1, 0, 1, 10, 0, 8, 4}, 
			{2, 0, 6, 10, 4, 26, 2}, 
			{3, 0, 4, 6, 2, 19, 1}, 
			{4, 0, 4, 8, 1, 3, 19}, 
			{5, 4, 8, 9, 19, 20, 28}, 
			{6, 5, 8, 9, 28, 23, 22},
			{7, 1, 5, 8, 23, 6, 5}, 
			{8, 1, 5, 7, 5, 21, 7},  
			{9, 1, 7, 10, 7, 8, 26}, 
			{10, 7, 10, 11, 26, 27, 29},  
			{11, 6, 10, 11, 29, 24, 25}, 
			{12, 2, 6, 11, 25, 13, 11},  
			{13, 2, 4, 6, 11, 18, 10}, 	
			{14, 2, 4, 9, 10, 20, 12}, 
			{15, 2, 3, 9, 12, 9, 16}, 
			{16, 3, 5, 9, 16, 23, 14}, 
			{17, 3, 5, 7, 14, 21, 15}, 
			{18, 3, 7, 11, 15, 27, 17}, 
			{19, 2, 3, 11, 17, 13, 9} 
		};
    return;
}



vector<vector<double>> base_spigoli;
vector<vector<double>> base_vertici;
    int new_id_v = mesh.M0D.size() + 1;
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




vector<vector<double>> vertici_originali = mesh.M0D; //vertici
vector<vector<double>> spigoli_originali = mesh.M1D;
vector<vector<double>> facce_originali = mesh.M2D; //facce

vector<vector<double>> tabella_confinamenti; //tab_conf
vector<vector<double>> vertici_duale; //v_duale
vector<vector<double>> spigoli_duale; //s_duale

void costruttore_tabella_confinamenti_3(vector<vector<double>>& F, vector<vector<double>>& T) {
    for (size_t i = 0; i < F.size(); ++i){
        T[i][0] = F[i][0]
        int primo_spigolo = F[i][4];
        int secondo_spigolo = F[i][5];
        int terzo_spigolo = F[i][6];
        int h = 1;
        for(const auto& j : F){
            if (F[j][4] == primo_spigolo || F[j][5] == primo_spigolo || F[j][6] == primo_spigolo ||
                F[j][4] == secondo_spigolo || F[j][5] == secondo_spigolo || F[j][6] == secondo_spigolo ||
                F[j][4] == terzo_spigolo || F[j][5] == terzo_spigolo || F[j][6] == terzo_spigolo) {
                    T[i][h] = F[j][0];
                    h++;
                }
        }
    }
    return;
}

void baricentro_3(vector<vector<double>>& V, vector<vector<double>>& F, vector<vector<double>>& VD){
    for (size_t i = 0; i < F.size(); ++i){
        int v1 = F[i][1]; //vertice 1
        int v2 = F[i][2]; //vertice 2
        int v3 = F[i][3]; //vertice 3                                           
        VD[i][0] = i + 1; //id vertici per il duale
        VD[i][1] = (V[v1 + 1][1] + V[v2 + 1][1] + V[v3 + 1][1])/3;
        VD[i][2] = (V[v1 + 1][2] + V[v2 + 1][2] + V[v3 + 1][2])/3;
        VD[i][3] = (V[v1 + 1][3] + V[v2 + 1][3] + V[v3 + 1][3])/3;
        VD[i][4] = F[i][0]; //a che faccia appartiene il baricentro
    }
    return;
}

void costruttore_duale_3(vector<vector<double>>& tab_conf, vector<vector<double>>& vertici, vector<vector<double>>& facce, vector<vector<double>>& v_duale, vector<vector<double>> s_duale){
    costruttore_tabella_confinamenti_3(facce, tab_conf);
    baricentro_3(vertici, facce, v_duale);
    int k = 1;
    for (size_t i = 0; i < v_duale.size(); ++i){
        int id_baricentro = v_duale[i][0];
        int id_faccia = v_duale[i].back();
        auto punta = find_if(tab_conf.begin(), tab_conf.end(), [id_faccia](const vector<int>& riga) {
            return !riga.empty() && riga[0] == id_faccia; //dentro il find_if c'è una cosiddetta funzione lambda (funzioni usa e getta), vedi teoria da qualche parte
        });
        //punta è un puntatore alla riga d'interesse.
        vector<double> riga_utile = *punta; //la riga relativa alla faccia del baricentro i-esimo nella matrice di prossimita (quella in cui ci sono gli id delle facce vicine)
        
        auto p_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[1]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[1]; // (= Primo Vertice Confinante)
        });
        int id_primo = (*p_v_c)[0];

        auto s_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[2]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[2]; // (= Secondo Vertice Confinante)
        });
        int id_secondo = (*s_v_c)[0];

        auto t_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[3]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[3]; // (= Terzo Vertice Confinante)
        });
        int id_terzo = (*t_v_c)[0];



        //ora inserisco i nuovi spigoli nella matrice spigoli del duale, solo dopo aver controllato che non ci siano già
        auto verifica = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_primo), max(id_baricentro,id_primo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_primo) && riga[2] == max(id_baricentro,id_primo);
        });
        
        if (verifica != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_primo), max(id_baricentro,id_primo)})
            k++;
        }

        auto verifica1 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_secondo), max(id_baricentro,id_secondo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_secondo) && riga[2] == max(id_baricentro,id_secondo);
        });
        
        if (verifica1 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_secondo), max(id_baricentro,id_secondo)})
            k++;
        }

        auto verifica2 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_terzo), max(id_baricentro,id_terzo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_terzo) && riga[2] == max(id_baricentro,id_terzo);
        });
        
        if (verifica2 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_terzo), max(id_baricentro,id_terzo)})
            k++;
        }
    }
}

void costruttore_tabella_confinamenti_4(vector<vector<double>>& F, vector<vector<double>>& T) {
    for (size_t i = 0; i < F.size(); ++i){
        T[i][0] = F[i][0]
        int primo_spigolo = F[i][5];
        int secondo_spigolo = F[i][6];
        int terzo_spigolo = F[i][7];
        int quarto_spigolo = F[i][8]
        int h = 1;
        for(const auto& j : F){
            if (F[j][5] == primo_spigolo || F[j][6] == primo_spigolo || F[j][7] == primo_spigolo ||
                F[j][5] == secondo_spigolo || F[j][6] == secondo_spigolo || F[j][7] == secondo_spigolo ||
                F[j][5] == terzo_spigolo || F[j][6] == terzo_spigolo || F[j][7] == terzo_spigolo ||
                F[j][5] == quarto_spigolo || F[j][6] == quarto_spigolo || F[j][7] == quarto_spigolo) {
                    T[i][h] = F[j][0];
                    h++;
                }
        }
    }
    return;
}

void baricentro_4(vector<vector<double>>& V, vector<vector<double>>& F, vector<vector<double>>& VD){
    sort(F.begin(), F.end(), [](const int& a, const int& b){
        return a[0] < b[0];
    }),
    for (size_t i = 0; i < F.size(); ++i){
        int v1 = F[i][1]; //vertice 1
        int v2 = F[i][2]; //vertice 2
        int v3 = F[i][3]; //vertice 3                                           occhio che le facce se create post triangolazione non sono ordinate!!
        int v4 = F[i][4]; //vertice 4
        VD[i][0] = i + 1; //id vertici per il duale

        VD[i][1] = (V[v1 + 1][1] + V[v2 + 1][1] + V[v3 + 1][1] + V[v4 + 1][1])/4;
        VD[i][2] = (V[v1 + 1][2] + V[v2 + 1][2] + V[v3 + 1][2] + V[v4 + 1][2])/4;
        VD[i][3] = (V[v1 + 1][3] + V[v2 + 1][3] + V[v3 + 1][3] + V[v4 + 1][3])/4;
        VD[i][4] = F[i][0]; //a che faccia appartiene il baricentro
    }
    return;
}

void costruttore_duale_4(vector<vector<double>>& tab_conf, vector<vector<double>>& vertici, vector<vector<double>>& facce, vector<vector<double>>& v_duale, vector<vector<double>> s_duale){
    costruttore_tabella_confinamenti_4(facce, tab_conf);
    baricentro_4(vertici, facce, v_duale);
    int k = 1;
    for (size_t i = 0; i < v_duale.size(); ++i){
        int id_baricentro = v_duale[i][0];
        int id_faccia = v_duale[i].back();
        auto punta = find_if(tab_conf.begin(), tab_conf.end(), [id_faccia](const vector<int>& riga) {
            return !riga.empty() && riga[0] == id_faccia; //dentro il find_if c'è una cosiddetta funzione lambda (funzioni usa e getta), vedi teoria da qualche parte
        });
        //punta è un puntatore alla riga d'interesse.
        vector<double> riga_utile = *punta; //la riga relativa alla faccia del baricentro i-esimo nella matrice di prossimita (quella in cui ci sono gli id delle facce vicine)
        
        auto p_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[1]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[1]; // (= Primo Vertice Confinante)
        });
        int id_primo = (*p_v_c)[0];

        auto s_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[2]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[2]; // (= Secondo Vertice Confinante)
        });
        int id_secondo = (*s_v_c)[0];

        auto t_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[3]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[3]; // (= Terzo Vertice Confinante)
        });
        int id_terzo = (*t_v_c)[0];

        auto q_v_c = find_if(v_duale.begin(), v_duale.end(), [riga_utile[4]](const vector<int>& riga) {
            return !riga.empty() && riga.back() == riga_utile[4]; // (= Terzo Vertice Confinante)
        });
        int id_quarto = (*q_v_c)[0];



        //ora inserisco i nuovi spigoli nella matrice spigoli del duale, solo dopo aver controllato che non ci siano già
        auto verifica1 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_primo), max(id_baricentro,id_primo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_primo) && riga[2] == max(id_baricentro,id_primo);
        });
        
        if (verifica1 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_primo), max(id_baricentro,id_primo)})
            k++;
        }


        auto verifica2 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_secondo), max(id_baricentro,id_secondo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_secondo) && riga[2] == max(id_baricentro,id_secondo);
        });
        
        if (verifica2 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_secondo), max(id_baricentro,id_secondo)})
            k++;
        }


        auto verifica3 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_terzo), max(id_baricentro,id_terzo)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_terzo) && riga[2] == max(id_baricentro,id_terzo);
        });
        
        if (verifica3 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_terzo), max(id_baricentro,id_terzo)})
            k++;
        }

        auto verifica4 = find_if(s_duale.begin(), s_duale.end(), [min(id_baricentro,id_quarto), max(id_baricentro,id_quarto)](const vector<int>& riga) {
            return riga.size() >= 3 && riga[1] == min(id_baricentro,id_quarto) && riga[2] == max(id_baricentro,id_quarto);
        });
        
        if (verifica4 != s_duale.end()){
            s_duale.push_back({k, min(id_baricentro,id_quarto), max(id_baricentro,id_quarto)})
            k++;
        }
    }
}