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

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void swap_val(int& a, int& b) {  //dati due valori fa lo swap
		int temp = a;
		a = b;
		b = temp;
	}
	
	vector<int> topological_Calc_I(const int q, const int b, const int c)
	{
		vector<int> risultato(3); // inizializza un vettore di dimensione 3 con tutti i valori pari a -1 
	
		int T = 0;
		T = b * b + b * c + c * c;
		int V = 0;
		int E = 0;
		int F = 0;
	
		if (q == 3) {
			V = 2 * T + 2;
			E = 6 * T;
			F = 4 * T;
		}
		else if (q == 4) {
			V = 4 * T + 2;
			E = 12 * T;
			F = 8 * T;
		}
		else {
			V = 10 * T + 2;
			E = 30 * T;
			F = 20 * T;
		}
	
		risultato[0] = V;  // il primo elemento corrisponde al numero di vertici V 
		risultato[1] = E;  // il secondo elemento corrisponde al numero di spigoli/lati E 
		risultato[2] = F;  // il terzo elemento corrisponde al numero di facce F 
	
		return risultato;  // Restituisce il vettore con i valori corrispondenti
	}
	
	//funzione che calcola i valori di V,E,F per la triangolazione di tipo 2
	vector<int> topological_Calc_II(const int b, const int q)
	{
		vector<int> output(3); // inizializza un vettore di dimensione 3 con tutti i valori pari a -1 
		int V = 0;
		int E = 0;
		int F = 0;
		
		if (q == 3) {
			V = 4 + 6 * (2 * b - 1) + static_cast<int>(round(4 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 6 * 2 * b + static_cast<int>(round(4 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 4 * (3 * b * b + 3 * b);
		}
		else if (q == 4) {
			V = 6 + 12 * (2 * b - 1) + static_cast<int>(round(8 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 12 * 2 * b + static_cast<int>(round(8 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 8 * (3 * b * b + 3 * b);
		}
		else {
			V = 12 + 30 * (2 * b - 1) + static_cast<int>(round(20 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 30 * 2 * b + static_cast<int>(round(20 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 20 * (3 * b * b + 3 * b);
		}
 
		//inserisco i valori calcolati all'interno del vettore output inizializzato all'inizio
		output[0] = V;  
		output[1] = E;  
		output[2] = F;
		
		return output; //restituisce il vettore calcolato
	}

	//funzione che calcola i valori di V,E,F del poliedro triangolato considerati i duplicati
	// dimension = vettore che contiene i valori di V,E,F del poliedro triangolato senza duplicati
	vector<int> duplicates_Calc(const int q, const int b, const int c, const vector<int>& dimension)
	{
		vector<int> valori_fin(3); 
		int sud_livello = 0; //variabile che suddivide il livello
		sud_livello = b + c;
		int V = dimension[0];
		int E = dimension[1];
		
		if (q == 3) {
			V += 2*4 + 6*(sud_livello - 1);
			E += 6*sud_livello;
		}
		else if (q == 4) {
			V += 3*6 + 12*(sud_livello - 1);
			E += 12*sud_livello;
		}
		else {
			V += 4*12 + 30*(sud_livello - 1);
			E += 30* sud_livello;
		}
		//inserisco i valori calcolati all'interno del vettore output inizializzato all'inizio
		valori_fin[0] = V;  
		valori_fin[1] = E;  
		valori_fin[2] = dimension[2]; //il numero delle facce non varia
		
		return valori_fin;
	}
	
}