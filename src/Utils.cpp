#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;

vector<vector<double>> vertici_tetraedro = {
    {  1, invSqrt3,  invSqrt3,  invSqrt3 },   
    {  2, invSqrt3, -invSqrt3, -invSqrt3 },   
    {  3, -invSqrt3,  invSqrt3, -invSqrt3 },   
    {  4, -invSqrt3, -invSqrt3,  invSqrt3 }
};


namespace PolygonalLibrary

{

bool importmesh(PolygonalMesh& mesh, vector<vector<double>>& A, vector<vector<int>>& B, vector<vector<int>>& C, vector<vector<int>>& D){ //correggi, incompleto
    unsigned int i=0;
    if(ImportCell0Ds(mesh, A)){i++}
    if(ImportCell1Ds(mesh, B)){i++}
    if(ImportCell2Ds(mesh, C, D)){i++}

    if(i!=3){return false}
    return true;
}

bool ImportCell0Ds(PolygonalMesh& mesh, vector<vector<double>>& A){
    mesh.M0D = A;
    return true;
}

bool Importcell1Ds(PolygonalMesh& mesh, vector<vector<int>>& A){
    mesh.M1D = A;
    return true;
}

bool ImportCell2Ds(PolygonalMesh& mesh, vector<vector<double>>& V, vector<vector<double>>& S){
    mesh.M2D_vertici = V;
    mesh.M2D_spigoli = S;
    return true;
}

}

vector<vector<double>> linspace(unsigned int indice_a, unsigned int indice_b, double x_a, double x_b, double y_a, double y_b, double z_a, double z_b, unsigned int N, unsigned int indice_min, unsigned int marker) {
    vector<vector<double>> result(N+1, vector<double>(3, 0.0));
    indice_min++;
    if (N < 0) return result;
    if (N = 0) {
        result[0].push_back(indice_a);
        result[0].push_back(x_a);
        result[0].push_back(y_a);
        result[0].push_back(z_a)
        return result;
    }
    result.front().push_back(indice_a, x_a, y_a);
    for (int i = 1; i < N; ++i){
        result[i].push_back(indice_min)
        indice_min++;
    }
    result.back().push_back(indice_b, x_b, y_b);
    double step_x = (x_b - x_a) / (N);

    for (int i = 0; i < N + 1; ++i) {
        result[i].push_back(x_a + i*step);
    }

    double step_y = (y_b - y_a) / (N);
    for (int i = 0; i < N + 1; ++i) {
        result[i].push_back(y_a + i * step);
    }

     double step_z = (z_b - z_a) / (N);
    for (int i = 0; i < N + 1; ++i) {
        result[i].push_back(z_a + i*step);
    }

    for (int i = 0; i < N + 1; i++){
        result[i].push_back(marker);
    }
    return result;
}

void controllo_ordine(vector<vector<double>>& matrice1_faccia, vector<vector<double>>& matrice2_faccia, vector<vector<double>>& matrice3_faccia){   
    for (unsigned int i : vettore1_faccia){
        if (matrice1_faccia[i][0] == matrice2_faccia[i][0]){
            return;
        }
    }
    reverse(matrice2_faccia.begin(), matrice2_faccia.end());
    return;
}

vector<double> creatore_vertici(unsigned int& i, double& x0, double& y0, double& z0, double& x1, double& y1, double& z1, double& x2, double& y2, double& z2){ 
    double diff_a = x1 - x0; 
    double diff_b = y1 - y0;
    double diff_c = z1 - z0;
    vector<double> nuovo_vertice = { i + 1, x2 + diff_a, y2 + diff_b, z2 + diff_c, 0 };
    i++;
    return nuovo_vertice;
}

void check_sovrapposizione(int& k, double& x, double& y, double& z, vector<vector<double>> matrice3){
        if (matrice3[k+1][1] ==  x && matrice3[+1][2] ==  y && matrice3[k+1][3] ==  z){
            return;
        }
    reverse(matrice2.begin(), matrice2.end());
    return;
}


bool mesh_tetraedro_I(PolygonalMesh& mesh, unsigned int b){
    vector<vector<double>> spigoli = mesh.M1D;
    vector<vector<double>> vertici = mesh.M0D;
    unsigned int dim = vertici.size();
    unsigned int l_righe = spigoli.size();
    unsigned int marker = 0;
    vector<vector<double>> nuovi_vertici;
    unsigned int max_indice = vertici.back().front(); //ricordati di gestire la conversione tipo di variabile
    for (int i = 0, i < l_righe, i++){
        unsigned int indice_partenza = spigoli[i][1];
        unsigned int indice_arrivo = spigoli[i][2];
        for (int i = 0, i < dim, i++){
            if (vertici[i][0] == indice_partenza){
                double x_partenza = vertici[i][1];
                double y_partenza = vertici[i][2];
                double z_partenza = vertici[i][3];
            }
        }
        for (int i = 0, i < dim, i++){
            if (vertici[i][0] == indice_arrivo){
                double x_arrivo = vertici[i][1];
                double y_arrivo = vertici[i][2];
                double z_arrivo = vertici[i][3];
            }
        }
        vector<vector<double>> divisione_lato1 = linspace(spigoli[i][1],spigoli[i][2], x_partenza, x_arrivo, y_partenza, y_arrivo, z_partenza, z_arrivo, b, max_indice, marker);
        nuovi_vertici.insert(nuovi_vertici.end(), divizione_lato.begin(), divisione_lato.end());
        marker++;        
    }

    vector<vector<double> vertici_duale; 
    
    vector<vector<double>> facce = mesh.M2D;
    vector<vector<double>> nuove_facce;
    vector<vector<double>> nuovi_spigoli;
    unsigned int n_facce = facce.size();
    for (int i = 0, i < n_facce, i++){
        unsigned int id_vertice1 = facce[i][1];
        unsigned int id_vertice2 = facce[i][2];
        unsigned int id_vertice3 = facce[i][3];
        unsigned int id_faccia1 = facce[i][4];
        unsigned int id_faccia2 = facce[i][5];
        unsigned int id_faccia3 = facce[i][6];
        vector<vertor<double>> = matrice1_faccia;
        vector<vertor<double>> = matrice2_faccia;
        vector<vertor<double>> = matrice3_faccia;
        for (auto i : nuovi_vertici){
            if (nuovi_vertici[i][4] == id_vertice1){
                matrice1_faccia.push_back(nuovi_vertici[i]);
            }
             if (nuovi_vertici[i][4] == id_vertice2){
                matrice2_faccia.push_back(nuovi_vertici[i]);
            }
             if (nuovi_vertici[i][4] == id_vertice3){
                matrice3_faccia.push_back(nuovi_vertici[i]);
            }
        }
        controllo_ordine(matrice1_faccia, matrice2_faccia, matrice3_faccia);

        int indice_vertice0 = matrice1_faccia[0][0];
        int indice_vertice1 = matrice1_faccia[1][0];
        int indice_vertice2 = matrice2_faccia[1][0];

        /*vector<vector<double> matrice_dati;
        
        for (unsigned int i = 0; i < 2; i++){
            for (unsigned int j = 1; j < 4; j++){
                matrice_dati.push_back(matrice1_faccia[i][j]);
            }
        }

        for (unsigned int i = 1; i < 4; i++){
            matrice_dati.push_back(matrice2_faccia[1][i]);
        }*/

        vector<double> punto_creato = {max_indice, x11 - x00 + x22, y11 - y00 + y22, z11 - z00 + z22};
        max_indice++;
        
        double x00 = matrice1_faccia[0][1];
        double y00 = matrice1_faccia[0][2];
        double z00 = matrice1_faccia[0][3];
        double x11 = matrice1_faccia[1][1];
        double y11 = matrice1_faccia[1][2];
        double z11 = matrice1_faccia[1][3];
        double x22 = matrice2_faccia[1][1];
        double y22 = matrice2_faccia[1][2];
        double z22 = matrice2_faccia[1][3];
        
        vector<double> punto_creato = {max_indice, x11 - x00 + x22, y11 - y00 + y22, z11 - z00 + z22};
        max_indice++;
        
        unsigned int index_facce = 0;
        vector<double> quarto_punto;
        vector<int> dati_faccia;
        vector<double> inserimento;
        vector<vector<double>> matrice_di_iterazione;
        vector<vector<double>> facce_provvisorie;
        int punto = matrice2_faccia[1][0];
        int b1 = 5;
        for (unsigned int j = 0, j < b - 1, j++){

            bool fine = false;
            unsigned int index = matrice2_faccia[j+1][0];

            for (unsigned int i = 0, i < b1 , i++){              
                nuovi_spigoli[3*i].push_back(3*i);
                nuovi_spigoli[3*i].push_back(matrice1_faccia[i][0]);
                nuovi_spigoli[3*i].push_back(matrice1_faccia[i+1][0]);  
                nuovi_spigoli[3*i+1].push_back(3*i + 1);
                nuovi_spigoli[3*i+1].push_back(matrice1_faccia[i][0]);
                nuovi_spigoli[3*i+1].push_back(index); 
                nuovi_spigoli[3*i+2].push_back(3*i + 2);
                nuovi_spigoli[3*i+2].push_back(matrice1_faccia[i+1][0]);
                nuovi_spigoli[3*i+2].push_back(index);                 

                dati_faccia = {index_facce, indice_vertice0, indice_vertice1, indice_vertice2, nuovi_spigoli[3*i][0], nuovi_spigoli[3*i + 1][0], nuovi_spigoli[3*i + 2][0]}
                nuove_facce.push_back(dati_faccia);
                dati_faccia.clear();
                index_facce++;

                if (i < b1 - 1){
                    quarto_punto = creatore_vertici(max_indice, x00 , y00, z00, x11, y11, z11, x22, y22, z22, 0);
                    if (i == b1 - 2){
                    quarto_punto[0] = matrice3_faccia[i+1][0];
                    fine = true;
                    }
                } 

                if (i == 1){
                    matrice_di_iterazione.push_back(matrice2_faccia[j+1]);
                } else {
                    matrice_di_iterazione.push_back(quarto_punto);
                }

                if (i > 0){
                    inserimento = {index_facce, matrice_di_iterazione[i][0], matrice1_faccia[i][0], matrice_di_iterazione[i-1][0], nuovi_spigoli[3*i-1][0], nuovi_spigoli[3*i+1][0]};
                    facce_provvisorie.push_back(inserimento);
                    inserimento.clear();
                    index_facce++;
                }
                
                if (fine == false){
                    dati_faccia = {index_facce, quarto_punto[0], indice_vertice1, indice_vertice};
                    index_facce++;
                }

                if (i == 0){
                    vector<double> punto_importante = quarto_punto;
                }
                
                if (i < b1 - 1){
                    indice_vertice0 = matrice1_faccia[i+1][0];
                    indice_vertice1 = matrice1_faccia[i+2][0];
                    indice_vertice2 = quarto_punto[0];
                    x00 = matrice1_faccia[i+1][1];
                    y00 = matrice1_faccia[i+1][2];
                    z00 = matrice1_faccia[i+1][3];
                    x11 = matrice1_faccia[i+2][1];
                    y11 = matrice1_faccia[i+2][2];
                    z11 = matrice1_faccia[i+2][3];
                    x22 = quarto_punto[1];
                    y22 = quarto_punto[2];
                    z22 = quarto_punto[3];
                } else {
                    indice_vertice0.clear();
                    indice_vertice1.clear();
                    indice_vertice2.clear();
                }
                
                index = quarto_punto[0];
                if (fine == false){
                    nuovi_vertici.push_back(quarto_punto);
                }

            }

            matrice1_faccia.clear();
            matrice1_faccia = matrice_di_iterazione;
            matrice_di_iterazione.clear();

            for (int k : facce_provvisorie){
                facce_provvisorie[k].push_back(matrice1_faccia[k]);
                nuove_facce.push_back(facce_provvisorie[k]);
            }

            facce_provvisorie.clear();

            x00 = matrice2_faccia[j][1];
            y00 = matrice2_faccia[j][2];
            z00 = matrice2_faccia[j][3];
            x11 = punto_importante[1];
            y11 = punto_importante[2];
            z11 = punto_importante[3];
            x22 = matrice2_faccia[j+1][1];
            y22 = matrice2_faccia[j+1][2];
            z22 = matrice2_faccia[j+1][3];
            
            b1 = b1 - 1;
        }
        nuovi_spigoli.push_back()

    }


}


bool duale(){
    
}
