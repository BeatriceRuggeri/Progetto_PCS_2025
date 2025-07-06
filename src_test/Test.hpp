#pragma once

#include <gtest/gtest.h>
#include <iostream>
#include "Eigen/Eigen"
#include <vector>

#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace Eigen;
using namespace PolyhedralLibrary



namespace PolyhedralTests
{



// .PHYSICALLY INPUT AS GROUND TRUTH MY DUDES.


	
    TEST(Polyhedral_Test,Test_Formula_VEF)
	{
	int b = 2;
	int c = 0;
	
	//p=3,q=3 tetrahedro
	int q1 = 3;
	vector<int> expected1 = {10, 24, 16}; //expected value
	vector<int> output1 = VEF_Calc(q1, b, c);
	EXPECT_EQ(expected1, output1);
	
	//p=3, q=4 Ottaedro
	int q2 = 4;
	vector<int> expected2 = {18, 48, 32};
	vector<int> output2 = VEF_Calc(q2, b, c);
	EXPECT_EQ(expected2, output2);
	
    //p=3, q=5 Icosaedro 
	int q3 = 5;
	vector<int> expected3 = {42, 120, 80};
	vector<int> output3 = VEF_Calc(q3, b, c);
	EXPECT_EQ(expected3, output3);
	}
	
	





   TEST(Polyhedral_Test, Test_Duplicates)
   {
	int b = 2;
	int c = 0;
	
	//p=3, q=3 ()
	int q1 = 3;
	vector<int> dim1 = {10, 24, 16};
	vector<int> expected1 = {24, 36, 16};
	vector<int> output1 = duplicates(q1, b, c, dim1);
	EXPECT_EQ(expected1, output1);
	
	//p=3, q=4
	int q2 = 4;
	vector<int> dim2 = {18, 48, 32};
	vector<int> expected2 = {48, 72, 32};
	vector<int> output2 = duplicates(q2, b, c, dim2);
	EXPECT_EQ(expected2, output2);
	
	//p=3, q=5
	int q3 = 5;
	vector<int> dim3 = {42, 120, 80};
	vector<int> expected3 = {120, 180 , 80};
	vector<int> output3 = duplicates(q3, b, c, dim3);
	EXPECT_EQ(expected3, output3);
   }

//TEST TRIANGOLAZIONE

  TEST(Polyhedral_Test, Test_Triangulation_Tetrahedron)
  {
	  
//definiamo le mesh

	PolyhedralMesh meshExpected;
	PolyhedralMesh meshOutput; //quello che carichiamo  dall'algoritmo 

//******************* vertici expected *************************************** Compatibilità (?** vector vector?)
	meshExpected.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

	meshExpected.Cell0DsCoordinates = MatrixXd::Zeros(3, 10); //questa l'abbiamo cambiata, prima non era inizzializata 
	
	//cioè i punti che vado a controllare (Punti medi di ogni lato che generiamo nella triangolazione)
	meshExpected.Cell0DsCoordinates.col(0)  <<  0, 0, 0.57735;
	meshExpected.Cell0DsCoordinates.col(1)  <<  -0.57735, -0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(2)  <<  0, -0.57735, 0;
	meshExpected.Cell0DsCoordinates.col(3)  <<  -0.57735, 0, 0;
	meshExpected.Cell0DsCoordinates.col(4)  <<  -0.57735,  0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(5)  <<  0, 0, -0.57735;
	meshExpected.Cell0DsCoordinates.col(6)  <<  0, 0.57735, 0;
	meshExpected.Cell0DsCoordinates.col(7)  <<  0.57735,  -0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(8)  <<  0.57735, 0, 0;
	meshExpected.Cell0DsCoordinates.col(9)  <<  0.57735, 0.57735, 0.57735;

	meshExpected.Cell0DsFlag = {};

//******************* lati e spigoli expected ***************************************
	meshExpected.Cell1DsId ={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

	meshExpected.Cell1DsExtrema = MatrixXi(24, 2);
	meshExpected.Cell1DsExtrema << 
		0, 6,
		0, 3,
		3, 6,
		0, 8,
		0, 9,
		2, 8,
		0, 2,
		0, 1,
		1, 2,
		2, 3,
		1, 3,
		2, 7,
		2, 5,
		3, 5,
		3, 4,
		4, 5,
		5, 6,
		4, 6,
		5, 7,
		7, 8,
		5, 8,
		6, 8,
		8, 9,
		6, 9;

	meshExpected.Cell1DsFlag = {};
	
	meshExpected.Cell1DsOriginalFlag = {1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0}; //boundary or not
	
//******************* faccie expected ***************************************
	meshExpected.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

	meshExpected.Cell2DsVertices = {
		{9, 0, 6},
		{0, 1, 3},
		{0, 3, 6},
		{6, 3, 4},
		{9, 8, 0},
		{8, 7, 2},
		{8, 2, 0},
		{0, 2, 1},
		{1, 2, 3},
		{2, 7, 5},
		{2, 5, 3},
		{3, 5, 4},
		{4, 5, 6},
		{5, 7, 8},
		{5, 8, 6},
		{6, 8, 9}
	};

	meshExpected.Cell2DsEdges = {
		{4, 0, 23},
		{7, 10, 1},
		{1, 2, 0},
		{2, 14, 17},
		{22, 3, 4},
		{19, 11, 5},
		{5, 6, 3},
		{6, 8, 7},
		{8, 9, 10},
		{11, 18, 12},
		{12, 13, 9},
		{13, 15, 14},
		{15, 16, 17},
		{18, 19, 20},
		{20, 21, 16},
		{21, 22, 23}
	};

//******************* polyhedra expected ***************************************

	meshExpected.Cell3DsId = {0}; 
	meshExpected.NumCells0Ds = {10};
	meshExpected.NumCells1Ds = {24};
	meshExpected.NumCells2Ds = {16};
	meshExpected.Cell3DsVertices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	meshExpected.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
	meshExpected.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
	
	
	
//******************* polyhedra calcolata dall'algoritmo ***************************************

    PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	//generateTetrahedron(mesh); //algo pietro
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	
	//manca fare queste funzioni
	
	vector<int> dim = topological_Calc_I(q, b, c);
	vector<int> dim_duplicated = duplicates_Calc(q, b, c, dim);
	
	
	
//QUESTO È DA ASSOLUTAMENTE CAMBIARE


/*

funzionepietro(mesh,meshTriangolazione,b,c,dim_duplicates); in mesh suddivide e triangola, salva in mesh triangolazione
pulire_edges
pulire_vertici
NuovaMesh(meshTriangolazione,MeshOutput,dim);



*/

	tri_build_I(mesh, meshTriangolazione, b, c, dim_duplicates); //Chiamo la funzione di triangolazione e utilizzo i q,b,c come input)
	v_duplicates_rm(meshTriangolazione); //pulisce duplicates
	l_duplicates_rm(meshTriangolazione); //pulisce vertici
	meshClean(meshTriangolazione, meshOutput, dim); //creo nuova mesh da meshTriangolazione in meshOuput
    Assemble_3D(meshOutput); //add geometric elements (points, faces, sides etc to a mesh or volume)-> once populated we can test it
 
 
 //******************* paragoniamo expected e output calcolato *************************************
 
//ID vert
    EXPECT_EQ(meshExpected.Cell0DsId, meshOutput.Cell0DsId); //ID matching test

//Coordinate vertici
    EXPECT_TRUE(meshExpected.Cell0DsCoordinates.isApprox(meshOutput.Cell0DsCoordinates, 1e-6)); //compares 3D position of vertices, floating point math is never exact so we use this

//ID spigoli
    EXPECT_EQ(meshExpected.Cell1DsId, meshOutput.Cell1DsId); //checks ID are the same

//Extrema spigoli
    EXPECT_TRUE(meshExpected.Cell1DsExtrema == meshOutput.Cell1DsExtrema); //Cell1DsExtrema : 2xN matrix storing start and end of vertices for each edge

//Vogliamo tenere questa?
    EXPECT_EQ(meshExpected.Cell1DsFlag, meshOutput.Cell1DsFlag); //comparing metadata flags: boundary markers and special attributes

//ID facce
    EXPECT_EQ(meshExpected.Cell2DsId, meshOutput.Cell2DsId); //face IDs must match exactly

//vertici faccie, we're assertng if our statements are equal or not, stops code if assertion fails, expect generates a non fatal failure.

    ASSERT_EQ(meshExpected.Cell2DsVertices.size(), meshOutput.Cell2DsVertices.size()); //same number of faces
    for (size_t i = 0; i < meshExpected.Cell2DsVertices.size(); ++i) { //now check that list of vertex indices matches for each face
        EXPECT_EQ(meshExpected.Cell2DsVertices[i], meshOutput.Cell2DsVertices[i]) << "Mismatch! Cell2DsVertices at face: " << i; //if any doesn't match then we get this message
    }

//spigoli faccie
    ASSERT_EQ(meshExpected.Cell2DsEdges.size(), meshOutput.Cell2DsEdges.size()); //same amount
    for (size_t i = 0; i < meshExpected.Cell2DsEdges.size(); ++i) { //scorre e controlla edges, bounds the face to be identical
        EXPECT_EQ(meshExpected.Cell2DsEdges[i], meshOutput.Cell2DsEdges[i]) << "Mismatch! Cell2DsEdges at face: " << i;
    }

//polyhedra
    EXPECT_EQ(meshExpected.Cell3DsId, meshOutput.Cell3DsId);
    EXPECT_EQ(meshExpected.Cell3DsVertices, meshOutput.Cell3DsVertices);
    EXPECT_EQ(meshExpected.Cell3DsEdges, meshOutput.Cell3DsEdges);
    EXPECT_EQ(meshExpected.Cell3DsFaces, meshOutput.Cell3DsFaces);
   } 
 
 //Quindi abbiamo controllato ID, connetività e approssimiamo coordinate, etc etc.
 
 
 }
 
 
//**************** Test Ordine Lati (fine di un arco coincida con quello successivo (e+1)%E


    TEST(Polyhedral_Test, Test_Ordine_Lati)
    { 
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	
	tetraedro_gen(mesh); 
	
	//define topological values
	
	int b = 2;
	int q = 3;
	int c = 0; 
	unsigned int maxFlag = numeric_limits<unsigned int>::max(); //to detect missing edges
	
	//calculate VEF values and how many duplicates to expect
	
	vector<int> dim=topological_Calc_I(q, b, c);
	vector<int> dim_duplicates = duplicates_Calc(q, b, c, dim); //maybe we don't need duplicates?
	
	//qua di nuovo abbiamo bisogno delle funzioni rispettive di triangolazione, che utilizziamo i dati controllati, generiamo, pulliamo a salviamo
	
	tri_build_I(mesh, meshTriangolazione, b, c, dim_duplicates);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	//Vertice E faccia == origine arco 

    //ciclo faccie mesh triangolazione, per controllare la consistenza geometrica: vertices.size()==edges.size()
     
    for (size_t f = 0; f < meshTriangolazione.Cell2DsId.size(); ++f) {
		const auto& edges = meshTriangolazione.Cell2DsEdges[f]; //lista dei lati di una faccia
        const auto& vertices = meshTriangolazione.Cell2DsVertices[f]; //lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici della faccia dopo scorrere
        
		
		ASSERT_EQ(vertices.size(), E) << "error: amount of vertices and sides don't coincide per face" << f; //which is a fatal assertion because it is an important consistency and avoid downstream bugs
		
//Edge Loop Consistency Check

		 for (size_t e = 0; e < E; ++e) {
//Here we've used the formulas indicated by the pdf file


/*

e = 0 → (0 + 1) % 4 = 1
e = 1 → 2 % 4 = 2
e = 2 → 3 % 4 = 3
e = 3 → 4 % 4 = 0  
//e qua torniamo, quindi scorriamo edges[1], edges[2], ... e torniamo a edges[0]

*/

			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            //this way we never exceed the range from 0 to E-1, % remainder/integer division
            
            //defensive initialization, each edge connects two vertices, edge: 42 connects v 10 and 12
            //da sopra: unsigned int maxFlag = numeric_limits<unsigned int>::max();
            //che mi fa: 4294967295   (if unsigned int is 32 bits) il massimo valore rappresentabile

            unsigned int currentEdgeOrigin = maxFlag; // COSÌ È CHIARAMENTE INVALIDA (non abbiamo ancora trovato i nostri endpoint)
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false; //utilizziamo questo bool e scorriamo
            for (unsigned int i = 0; i < meshTriangolazione.Cell1DsId.size(); i++) {
                if (currentEdge == meshTriangolazione.Cell1DsId[i]) {// we find
                    currentEdgeOrigin = meshTriangolazione.Cell1DsExtrema(i, 0); //we replace
                    currentEdgeEnd = meshTriangolazione.Cell1DsExtrema(i, 1);
                    foundCurrent = true; //change the flag
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            
            bool foundNext = false;
            for (unsigned int i = 0; i < meshTriangolazione.Cell1DsId.size(); i++) {
                if (nextEdge == meshTriangolazione.Cell1DsId[i]) {
                    nextEdgeOrigin = meshTriangolazione.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshTriangolazione.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext) << "Master Edge ID " << nextEdge << " not found in Cell1DsId!";

            // recupero il vertice della faccia
            unsigned int faceVertex = vertices[e]; 

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);

            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);

		}		
	}
}



//************** test area non nulla ******************

    TEST(Polyhedral_Test, Test_Area)
    {
	double eps = numeric_limits<double>::epsilon(); //smallest representable positive difference in double precision (aprox 2.26e-16)

	// mesh dalla triangolazione 
	PolyhedralMesh meshTriangolazione; //generiamo mesh
	PolyhedralMesh mesh;
	tetraedro_gen(mesh); //generiamo il nostro tetrahedro dentro mesh
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dim=topological_Calc_I(q, b, c);
	vector<int> dim_duplicates = duplicates_Calc(q, b, c, dim);
	
	//our usual process to be defined
	
	tri_build_I(mesh, meshTriangolazione, b, c, dim_duplicates);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	
	//ciclo triangoli
	for (size_t i = 0; i < meshTriangolazione.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangolazione.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); //controllo quantità vertici, cioè controllando se ha 3 vertici controllo che sia un triangolo effettivamente 
		
		// accedo alle coordinate dei vertici
		
		/*
		
		A = [x0,y0,z0]
        B = [x1,y1,z1]
        C = [x2,y2,z2]

		*/
		
		Vector3d A = meshTriangolazione.Cell0DsCoordinates.col(tri[0]); //three coord vector: A=[x0,y0,z0];
        Vector3d B = meshTriangolazione.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangolazione.Cell0DsCoordinates.col(tri[2]);
		
	//area calc, speriamo che sia maggiore di una tolleranza
		double area = 0.5 * ((B - A).cross(C - A)).norm(); //vector product
		EXPECT_GT(area, eps) << "Error: Null area (sotto tolleranza)" << i;
	  }	
    }
// se non funziona, output di Google Test:
//Error: Null area (sotto tolleranza) i=17, così veddiamo precisamente quale triangolo non è valido

    
//così avviamo compiutato l'area senza assumere che siano piani


//************** test lati non nulli *********************

    TEST(Polyhedral_Test, Test_Lati)
    {
	double eps = numeric_limits<double>::epsilon();

	// usual procedure, create two meshes and generate in mesh our tetrahedron
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	tetraedro_gen(mesh); // !!!!!!!!!!!
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	//calculate topological values
	vector<int> dim=VEF_Calc(q, b, c);
	vector<int> dim_duplicates = duplicates(q, b, c, dim);
	
	//the thing to really replace
	
	tri_build_I(mesh, meshTriangolazione, b, c, dim_duplicates);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	
	//iterazione lati mesh triangolazione
	 
	for (unsigned int i = 0; i < meshTriangolazione.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vec_start = meshTriangolazione.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vec_end = meshTriangolazione.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// le prendo
		Vector3d start_point = meshTriangolazione.Cell0DsCoordinates.col(vec_start);
        Vector3d end_point = meshTriangolazione.Cell0DsCoordinates.col(vec_end);
		
		// calcolo la norma della lunghezza del lato
		double length = (end_point - start_point).norm(); //euclidian length
		
		//controllo con google test
		
		EXPECT_GT(length, eps) << "error: null length (under tolerance)" << i;
		
	 }

	    
    }

//******** test lati duale non nulli *******************

TEST(Polyhedral_Test, Test_Lati_Duale){
	
	double eps = numeric_limits<double>::epsilon();
	PolyhedralMesh meshDual;
	PolyhedralMesh mesh;
	tetraedro_gen(mesh); // Prendiamo nuovamente il tetrahedro e lo mettiamo in mesh, per poi fare la mesh dual e salvarlo in meshDual
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	triangolazione_dual(q, b, c, mesh, meshDual); // dove abbiamo fatto questo?
	
	
	//iterazione su lati triangolate
	for (unsigned int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		//dove prendiamo colonna 0 e 1 per ogni i-esima riga
		int vec_start = meshDual.Cell1DsExtrema(i, 0); // indice del vertice di partenza
        int vec_end = meshDual.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// le metto in questo vector
		Vector3d start_point = meshDual.Cell0DsCoordinates.col(vec_start);
        Vector3d end_point = meshDual.Cell0DsCoordinates.col(vec_end);
		
		// calcolo la norma della lunghezza del lato
		double length = (end_point - start_point).norm();
		
		EXPECT_GT(length, eps) << "error: side with null length (under tolerance)" << i;
		
	}
}

//********* test lati triangolazione due non nulla


TEST(Polyhedral_Test, Test_Lati_Tri_Due){
	
	double eps = numeric_limits<double>::epsilon();
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	
	//Dual of a tetrahedron is also a tetrahedron though guys (it is a self dual quantity) 
	
	tetraedri_gen(mesh);
	//quindi se coincidono allora va bene l'algoritmo, e il test anche.
	
	int q = 3;
	int b = 2;
	
   triangolazione_II(q, b, mesh, meshTriangolazione); 
	
	
	//iterazione su triangolata
	for (unsigned int i = 0; i < meshTriangolazione.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vec_start = meshTriangolazione.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vec_end = meshTriangolazione.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d start_point = meshTriangolazione.Cell0DsCoordinates.col(vec_start);
        Vector3d end_point = meshTriangolazione.Cell0DsCoordinates.col(vec_end);
		
		// calcolo la norma della lunghezza del lato
		double length = (end_point - start_point).norm();
		
		
		//controlliamo che sia non degerene
		EXPECT_GT(length, eps) << "error: side with null length (under tolerance)" << i;
		
	}
}

//this way we have no collapsed edges in our mesh


//**** test triangolazione due area non nulla

TEST(Polyhedral_Test, Test_Area_Tri_Due)
{
	double eps = numeric_limits<double>::epsilon();

	PolyhedralMesh mesh; //define mesh
	tetraedro_gen(mesh); //generate tetrahedron inside the mesh
	
	PolyhedralMesh meshTriangolazione; //define a couple more meshes (intermediate meshes)
	PolyhedralMesh meshOutput; //our final mesh
	PolyhedralMesh meshTriangolazione_Due; //second triangolation mesh
	
	int q = 3;
	int b = 2;
	
	//calculate topologies
	
	vector<int> dim= topological_Calc_I(q, b, 0);
	vector<int> dim_duplicates = duplicates_Calc(q, b, 0, dim);
	
	//************************* !!!!!!!!!! our triangulation pipeline
	tri_build_I(mesh, meshTriangolazione, b, 0,  dim_duplicates);
	v_duplicates_rm(meshTriangolazione;
	l_duplicates_rm(meshTriangolazione);
	meshClean(meshTriangolazione, meshOutput, dim);
	
	vector<int> dim2 = topological_Calc_II(b, q); //we get the expected count of elements in our second triangulation
	
	
	/*
	
	#include <map>
#include <string>

std::map<std::string, int> ages;
ages["Alice"] = 30;
ages["Bob"] = 25;

if (ages.find("Alice") != ages.end()) {
    int a = ages["Alice"]; // = 30
}

	
	*/
	
	//std::map<Key, Value> so each key is unique to one value (balanced binary tree)
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> LF_map = buildMap_LF(meshOutput); //look up table mapping the edges to the faces they belong to.
	//edge (2,5) => [face 10, face 22] basically

	tri_build_II(meshOutput, meshTriangolazione, dim2, edgeToFacesMap); //our other function
	
	// ciclo su tutti i triangoli
	for (size_t i = 0; i < meshTriangolazione_Due.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangolazione_Due.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); // controllo se il triangolo ha 3 vertici
		
		// accedo alle coordinate dei vertici
		Vector3d A = meshTriangolazione_Due.Cell0DsCoordinates.col(tri[0]); //prendo la colonna che contiene le coordinate
        Vector3d B = meshTriangolazione_Due.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangolazione_Due.Cell0DsCoordinates.col(tri[2]);
		
		//area calc
		double area = 0.5 * ((B - A).cross(C - A)).norm(); // prodotto vettoriale
		EXPECT_GT(area, eps) << "Triangolo con area nulla o quasi nulla al triangolo " << i;
	}	
	
}


//********* test duale, da rifare pratticamente

/*we will control that:

-topology is correct
-dual vertices are the baricenters of faces
-edges are correctly connected to adjacent faces
-Faces are consistently oriented

*/

TEST(Polyhedral_Test, Test_Duale){
	//Data una mesh triangolata di un tetraedro, costruisco la sua mesh polyhedral duale e controllo topologia, baricentri, lati, orientazione faccie
	// mesh ottenuta utilizzando la funzione di triangolazioneMore actions
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	PolyhedralMesh meshOutput; //i.e: mesh output prima
	PolyhedralMesh meshDual;
	tetraedro_gen(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = topological_Calc_I(q, b, c);
	vector<int> dimensionDuplicated = duplicates_Calc(q, b, c, dimension);
	
	
	tri_build_I(mesh, meshTriangulated, b, c, dimensionDuplicated);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	meshClean(meshTriangolazione, meshOutout, dimension);//so then we put it in mesh final
	//now we use a map to pair a key 
	map <pair<unsigned int, unsigned int>, vector<unsigned int>> LF_map = buildMap_LF(meshFinal); //so for every edge v0,v1 we will record which face shares it
	
	topological_Calc_II(meshOutput, meshDual, LF_Map);
	
	//define tolerances, which we will use for floating point comparisons
	double eps = numeric_limits<double>::epsilon();
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	
//these will be our expected counts
	size_t e_vert_dual= meshFinal.Cell2DsId.size(); // 16 facce triangolate → 16 vertici nel duale
    size_t e_lati_dual= meshFinal.Cell1DsId.size(); // Calcolato dinamicamente, può essere 24 per subdivisionLevel=2
   
    // numero di vertici (id)
	// numero di lati(id)
   
   //so now we check our counts:
    EXPECT_EQ(meshDual.Cell0DsId.size(), e_vert_dual);
    EXPECT_EQ(meshDual.Cell1DsId.size(), e_lati_dual);
	
	// ogni vertice del poliedro originale genera una faccia nel duale, i.e face cardinality
	EXPECT_GE(meshDual.Cell2DsId.size(), 4);  // Minimo 4 se è un tetraedro chiuso
	
	// verifico se ogni faccia duale ha almento 3 vertici
	for (const auto& dualFace : meshDual.Cell2DsVertices) {
    	EXPECT_GE(dualFace.size(), 3);
    }
	
	// we gotta check that they don't connect to itself:
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		int a = meshDual.Cell1DsExtrema(i, 0); //prendo uno
		int b = meshDual.Cell1DsExtrema(i, 1); //prendo l'altro
		EXPECT_NE(a, b); // Mi verifica che a e b siano diversi 
    }
	
	// verifico che il calcolo del baricentro sia corretto
	// Per ogni faccia della mesh triangolata
    for (size_t faceId = 0; faceId < meshOutput.Cell2DsId.size(); ++faceId) {
        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
		
        // accedo alla lista dei vertici che compongono la faccia
        const vector<unsigned int>& faceVertices = meshOutput.Cell2DsVertices[faceId];

        // sommo le coordinate dei vertici della faccia
        for (unsigned int v_id : faceVertices) {
			sumX += meshOutput.Cell0DsCoordinates(0, v_id); // coordinata x del vertice v_idMore actions
            sumY += meshOutput.Cell0DsCoordinates(1, v_id); // coordinata y del vertice v_id
            sumZ += meshOutput.Cell0DsCoordinates(2, v_id); // coordinata z del vertice v_id
        }

        // calcolo il baricentro
        double n = (double) faceVertices.size(); // perché faceVertices.size() restituisce un tipo size_t quindi evito una conversione implicita
        double baryX = sumX / n;
        double baryY = sumY / n;
        double baryZ = sumZ / n;

        // estraggo la coordinata del vertice duale corrispondente (baricentro)
        double X = meshDual.Cell0DsCoordinates(0, faceId);
        double Y = meshDual.Cell0DsCoordinates(1, faceId);
        double Z = meshDual.Cell0DsCoordinates(2, faceId);

        EXPECT_NEAR(baryX, X, eps); //e finalmente paragoniamo ragaa
        EXPECT_NEAR(baryY, Y, eps);
        EXPECT_NEAR(baryZ, Z, eps);
    }

    // verifico che i lati siano costruiti nel modo corretto ---> devono connettere i baricentri adjacent to each other
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		Vector2i edge = meshDual.Cell1DsExtrema.row(i);
		// id delle facce del poliedro originale (vertici duali sono baricentri di facce originali)
        int f1 = edge(0); 
        int f2 = edge(1);
		 // controlla che le due facce originali condividano un edge
    bool foundCommonEdge = false;
    for (int e1 : meshFinal.Cell2DsEdges[f1]) {
        for (int e2 : meshFinal.Cell2DsEdges[f2]) {
            if (e1 == e2) {
                foundCommonEdge = true;
                break;
            }
        }
        if (foundCommonEdge) break;
    }
    EXPECT_TRUE(foundCommonEdge);
	}
	
	// verifico che i lati e le facce siano ordinate correttamente (Consistency check)
	//Tutti lati connessi consecutivamente, topologicamente validi
	 
	for (size_t f = 0; f < meshDual.Cell2DsId.size(); ++f) {
		const auto& edges = meshDual.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshDual.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
//We can't assume the edge is stored with their vertex, we have to find it and confirm it exists, ensure the loop is closed by checking that consecutive edges share a vertex and check that the face vertex list is CONSISTENT with edge connnectivity
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag; //similar procedure that we did above
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (currentEdge == meshDual.Cell1DsId[i]) {
                    currentEdge_start = meshDual.Cell1DsExtrema(i, 0);
                    currentEdge_end = meshDual.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
            
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (nextEdge == meshDual.Cell1DsId[i]) {
                    nextEdge_start = meshDual.Cell1DsExtrema(i, 0);
                    nextEdge_end = meshDual.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext);

            // Recupera il vertice della faccia
            unsigned int faceVertex = vertices[e]; // ID master del vertice della faccia

            bool condition = (currentEdge_start == nextEdge_start) ||
                             (currentEdge_start == nextEdge_end) ||
                             (currentEdge_end == nextEdge_start) ||
                             (currentEdge_end == nextEdge_end);
            EXPECT_TRUE(condition);
              
            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            bool condition_II = (face_vert == currentEdge_start) ||
            					(face_vert == currentEdge_end);
            EXPECT_TRUE(condition_II);
		}		
	}	 
}


//*********** test triangolazione due



TEST(Test_Polyhedral, Test_Lati_II)
{ 
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::tetraedron_gen(mesh);
	
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh meshOutput;
	PolyhedralMesh meshTriangolazione_II;
	
	int q = 3;
	int b = 2;
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dim = topological_Calc_I(q, b, 0);
	vector<int> dim_dupli = duplicates_Calc(q, b, 0, dime);
	tri_build_I(mesh, meshTriangolazione, b, 0,  dim_dupli);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	meshClean(meshTriangolazione, meshOutput, dim);
	
	vector<int> dim_II = topological_Calc_II(b, q);
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> LF_map = buildMap_LF(meshOutput);
	tri_build_II(meshOutput, meshTriangolazione_II, dim_II, LF_map);
	
	// per ogni faccia i lati sono ordinati in modo che la fine dell'arco e coincida con l'inizio dell'arco successivo (e+1)%E
	// il vertice e della faccia deve corrispondere all'origine dell'arco e

    // ciclo su tutte le facce della mesh triangolata 
    for (size_t f = 0; f < meshOutput.Cell2DsId.size(); ++f) {
		const auto& edges = meshOutput.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshOutput.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
        
		
		ASSERT_EQ(vertices.size(), E) << "Numero di vertici e di lati non corrispondono per faccia " << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdge_start = maxFlag;
			unsigned int currentEdge_end = maxFlag;
			unsigned int nextEdge_start = maxFlag;
			unsigned int nextEdge_end = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshOutput.Cell1DsId.size(); i++) {
                if (currentEdge == meshOutput.Cell1DsId[i]) {
                    currentEdge_start = meshOutput.Cell1DsExtrema(i, 0);
                    currentEdge_end = meshOutput.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshOutput.Cell1DsId.size(); i++) {
                if (nextEdge == meshOutput.Cell1DsId[i]) {
                    nextEdge_start = meshOutput.Cell1DsExtrema(i, 0);
                    nextEdge_end = meshOutput.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext) << "Master Edge ID " << nextEdge << " not found in Cell1DsId";

            // recupero il vertice della faccia
            unsigned int faceVertex = vertices[e]; 

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);

            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);

		}		
	}
}

  
  
  
  
//******************* test cammini minimi (shortest path) ***************************************
//******adapt*******

TEST(Polyhedral_Test, ShortestPath){
	
	PolyhedralMesh mesh;
	PolyhedralMesh meshTriangolazione; //fissata expected
	PolyhedralMesh meshOutput; 
	
	//chiamo generazione tetrahedron
	
	tetraedro_gen(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	int startVertexId = 0;
	int endVertexId = 7;
	
	vector<int> dim = topological_Calc_I(q, b, c);
	vector<int> dim_duplicates = duplicates_Calc(q, b, c, dim);
	
//************************

	tri_build_I(mesh, meshTriangolazione, b, c, dim_duplicates);
	v_duplicates_rm(meshTriangolazione);
	l_duplicates_rm(meshTriangolazione);
	meshClean(meshTriangolazione, meshOutput, dim);
	
	minimo_output expected(2, 1.63299); //these values!! ?
	MatrixXi adjMatrix = adjMatrix(meshOutput);
	minimo_output res = minimo_F_Dijkstra(meshOutput, adjMatrix, startVertexId, endVertexId);
	
	EXPECT_NEAR(res.num_lati, expected.num_lati, 1e-5) 
		<< "error, not the shortest path"
		<< " Expected: " << expected.num_lati << ", Final: " << res.num_lati;
		
	EXPECT_NEAR(res.length, expected.length, 1e-5) 
        << "error, not the shortest path"
        << " Expected: " << expected.length << ", Final: " << res.length;
	 
	}
 
}
