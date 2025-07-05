#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;


//???????????????????????

vector<vector<double>> vertici_originali = mesh.M0D; //vertici
vector<vector<double>> spigoli_originali = mesh.M1D;
vector<vector<double>> facce_originali = mesh.M2D; //facce

vector<vector<double>> tabella_confinamenti; //tab_conf
vector<vector<double>> vertici_duale; //v_duale
vector<vector<double>> spigoli_duale; //s_duale

void costruttore_tabella_confinamenti(vector<vector<double>>& F, vector<vector<double>>& T) {
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

void baricentro(vector<vector<double>>& V, vector<vector<double>>& F, vector<vector<double>>& VD){
    for (size_t i = 0; i < F.size(); ++i){
        int v1 = F[i][1]; //vertice 1
        int v2 = F[i][2]; //vertice 2
        int v3 = F[i][3]; //vertice 3                                           occhio che le facce se create post triangolazione non sono ordinate!!
        VD[i][0] = i + 1; //id vertici per il duale
        VD[i][1] = (V[v1 + 1][1] + V[v2 + 1][1] + V[v3 + 1][1])/3;
        VD[i][2] = (V[v1 + 1][2] + V[v2 + 1][2] + V[v3 + 1][2])/3;
        VD[i][3] = (V[v1 + 1][3] + V[v2 + 1][3] + V[v3 + 1][3])/3;
        VD[i][4] = F[i][0]; //a che faccia appartiene il baricentro
    }
    return;
}

void costruttore_duale(vector<vector<double>>& tab_conf, vector<vector<double>>& vertici, vector<vector<double>>& facce, vector<vector<double>>& v_duale, vector<vector<double>> s_duale){
    costruttore_tabella_confinamenti(facce, tab_conf);
    baricentro(vertici, facce, v_duale);
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

void costruttore_tabella_confinamenti(vector<vector<double>>& F, vector<vector<double>>& T) {
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

void baricentro(vector<vector<double>>& V, vector<vector<double>>& F, vector<vector<double>>& VD){
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

void costruttore_duale2(vector<vector<double>>& tab_conf, vector<vector<double>>& vertici, vector<vector<double>>& facce, vector<vector<double>>& v_duale, vector<vector<double>> s_duale){
    costruttore_tabella_confinamenti(facce, tab_conf);
    baricentro(vertici, facce, v_duale);
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