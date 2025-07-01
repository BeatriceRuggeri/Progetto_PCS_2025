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



	//Test formula VEF: uses VEF_Calc function
	
    TEST(Polyhedral_Test,Test_Formula_VEF)
	{
	int b = 2;
	int c = 0;
	
	//p=3,q=3 tetrahedro
	int q1 = 3;
	vector<int> expected1 = {10, 24, 16};
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
	
	



//TEST DUPLICATES: uses duplicates function

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
	PolyhedralMesh meshOutput; //quello che buta fuori l'algoritmo 

//******************* vertici expected *************************************** Compatibilità (?** vector vector?)
	meshExpected.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

	meshExpected.Cell0DsCoordinates = MatrixXd::Zeros(3, 10);
	
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
	
	vector<int> dim = VEF_Calc(q, b, c);
	vector<int> dim_duplicated = duplicates(q, b, c, dim);
	triangulateAndStore(mesh, meshTriangolazione, b, c, dim_duplicates);
	RemoveDuplicatedEdges(meshTriangolazione);
	RemoveDuplicatedVertices(meshTriangolazione);
	NewMesh(meshTriangolazione, meshOutput, dim);
    PopulateCell3D(meshOutput);
 
 
 //******************* paragoniamo expected e output calcolato *************************************
 
//ID vert
    EXPECT_EQ(meshExpected.Cell0DsId, meshOutput.Cell0DsId);

//Coordinate vertici
    EXPECT_TRUE(meshExpected.Cell0DsCoordinates.isApprox(meshOutput.Cell0DsCoordinates, 1e-6));

//ID spigoli
    EXPECT_EQ(meshExpected.Cell1DsId, meshOutput.Cell1DsId);

//Extrema spigoli
    EXPECT_TRUE(meshExpected.Cell1DsExtrema == meshOutput.Cell1DsExtrema);

//Flag spigoli
    EXPECT_EQ(meshExpected.Cell1DsFlag, meshOutput.Cell1DsFlag);

//ID facce
    EXPECT_EQ(meshExpected.Cell2DsId, meshOutput.Cell2DsId);

//vertici faccie
    ASSERT_EQ(meshExpected.Cell2DsVertices.size(), meshOutput.Cell2DsVertices.size());
    for (size_t i = 0; i < meshExpected.Cell2DsVertices.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsVertices[i], meshOutput.Cell2DsVertices[i]) << "Mismatch! Cell2DsVertices at face: " << i;
    }

//spigoli faccie
    ASSERT_EQ(meshExpected.Cell2DsEdges.size(), meshOutput.Cell2DsEdges.size());
    for (size_t i = 0; i < meshExpected.Cell2DsEdges.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsEdges[i], meshOutput.Cell2DsEdges[i]) << "Mismatch! Cell2DsEdges at face: " << i;
    }

//polyhedra
    EXPECT_EQ(meshExpected.Cell3DsId, meshOutput.Cell3DsId);
    EXPECT_EQ(meshExpected.Cell3DsVertices, meshOutput.Cell3DsVertices);
    EXPECT_EQ(meshExpected.Cell3DsEdges, meshOutput.Cell3DsEdges);
    EXPECT_EQ(meshExpected.Cell3DsFaces, meshOutput.Cell3DsFaces);
   } 
 
 
 
 }
 
 
//**************** Test Ordine Lati (fine di un arco coincida con quello successivo (e+1)%E


    TEST(Polyhedral_Test, Test_Ordine_Lati)
    { 
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	generateTetrahedron(mesh); //!!!!!!
	
	int b = 2;
	int q = 3;
	int c = 0; 
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dim=VEF_Calc(q, b, c);
	vector<int> dim_duplicates = duplicates(q, b, c, dim);
	triangulateAndStore(mesh, meshTriangolazione, b, c, dim_duplicates);
	RemoveDuplicatedVertices(meshTriangolazione);
	RemoveDuplicatedEdges(meshTriangolazione);
	//Vertice E faccia == origine arco 

    //ciclo faccie mesh triangolazione
    for (size_t f = 0; f < meshTriangolazione.Cell2DsId.size(); ++f) {
		const auto& edges = meshTriangolazione.Cell2DsEdges[f]; //lista dei lati di una faccia
        const auto& vertices = meshTriangolazione.Cell2DsVertices[f]; //lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici della faccia
        
		
		ASSERT_EQ(vertices.size(), E) << "error: amount of vertices and sides don't coincide per face" << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshTriangolazione.Cell1DsId.size(); i++) {
                if (currentEdge == meshTriangolazione.Cell1DsId[i]) {
                    currentEdgeOrigin = meshTriangolazione.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshTriangolazione.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
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
	double eps = numeric_limits<double>::epsilon();

	// mesh dalla triangolazione 
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	generateTetrahedron(mesh); // !!!!!!!!!!!!
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dim=VEF_Calc(q, b, c);
	vector<int> dim_duplicates = duplicates(q, b, c, dim);
	triangulateAndStore(mesh, meshTriangolazione, b, c, dim_duplicates);
	RemoveDuplicatedEdges(meshTriangolazione);
	RemoveDuplicatedVertices(meshTriangolazione);
	
	//ciclo triangoli
	for (size_t i = 0; i < meshTriangolazione.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangolazione.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); //controllo vertici
		
		// accedo alle coordinate dei vertici
		Vector3d A = meshTriangolazione.Cell0DsCoordinates.col(tri[0]); 
        Vector3d B = meshTriangolazione.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangolazione.Cell0DsCoordinates.col(tri[2]);
		
	//area calc
		double area = 0.5 * ((B - A).cross(C - A)).norm(); //vector product
		EXPECT_GT(area, eps) << "Error: Null area (sotto tolleranza)" << i;
	  }	
    }


//************** test lati non nulli *********************

    TEST(Polyhedral_Test, Test_Lati)
    {
	double eps = numeric_limits<double>::epsilon();

	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh mesh;
	generateTetrahedron(mesh); // !!!!!!!!!!!
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dim=VEF_Calc(q, b, c);
	vector<int> dim_duplicates = duplicates(q, b, c, dim);
	triangulateAndStore(mesh, meshTriangolazione, b, c, dim_duplicates);
	RemoveDuplicatedEdges(meshTriangolazione);
	RemoveDuplicatedVertices(meshTriangolazione);
	
	//iterazione lati mesh triangolazione 
	for (unsigned int i = 0; i < meshTriangolazione.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vec_start = meshTriangolazione.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vec_end = meshTriangolazione.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d start_point = meshTriangolazione.Cell0DsCoordinates.col(vec_start);
        Vector3d end_point = meshTriangolazione.Cell0DsCoordinates.col(vec_end);
		
		// calcolo la norma della lunghezza del lato
		double length = (end_point - start_point).norm();
		
		EXPECT_GT(length, eps) << "error: null length (under tolerance)" << i;
		
	 }

	    
    }

//******** test lati duale non nulli *******************

TEST(Polyhedral_Test, Test_Lati_Duale){
	
	double eps = numeric_limits<double>::epsilon();
	PolyhedralMesh meshDual;
	PolyhedralMesh mesh;
	generateTetrahedron(mesh); // !!!!!!!
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	TriangulationDual(q, b, c, mesh, meshDual); // dove abbiamo fatto questo?
	
	
	//iterazione su lati triangolate
	for (unsigned int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		int vec_start = meshDual.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vec_end = meshDual.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
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
	
	generateTetrahedron(mesh); //!!!!!!!
	
	int q = 3;
	int b = 2;
	
	Triangulation2(q, b, mesh, meshTriangolazione); //!!!!!!!!
	
	
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
		
		EXPECT_GT(length, eps) << "error: side with null length (under tolerance)" << i;
		
	}
}



//**** test triangolazione due area non nulla

TEST(Polyhedral_Test, Test_Area_Tri_Due)
{
	double eps = numeric_limits<double>::epsilon();

	PolyhedralMesh mesh;
	generateTetrahedron(mesh);
	
	PolyhedralMesh meshTriangolazione;
	PolyhedralMesh meshOutput;
	PolyhedralMesh meshTriangolazione_Due;
	
	int q = 3;
	int b = 2;
	
	vector<int> dim= VEF_Calc(q, b, 0);
	vector<int> dim_duplicates = duplicates(q, b, 0, dim);
	
	//************************* !!!!!!!!!!
	triangulateAndStore(mesh, meshTriangolazione, b, 0,  dim_duplicates);
	RemoveDuplicatedVertices(meshTriangolazione;
	RemoveDuplicatedEdges(meshTriangolazione);
	NewMesh(meshTriangolazione, meshOutput, dim);
	
	//************************ !!!!!!!!!
	
	vector<int> dim2 = CalculateDimension2(b, q);
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshOutput);
	triangulateAndStore2(meshOutput, meshTriangolazione, dim2, edgeToFacesMap);
	
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

TEST(Polyhedral_Test, Test_Duale){
	
	// mesh ottenuta utilizzando la funzione di triangolazioneMore actions
	PolyhedralMesh meshTriangulated;

	PolyhedralMesh mesh;
	PolyhedralMesh meshFinal;
	PolyhedralMesh meshDual;
	generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension = ComputePolyhedronVEF(q, b, c);
	vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RemoveDuplicatedEdges(meshTriangulated);
	RemoveDuplicatedVertices(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	map <pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
	
	CalculateDual(meshFinal, meshDual, edgeToFacesMap);
	
	double eps = numeric_limits<double>::epsilon();
	unsigned int maxFlag = numeric_limits<unsigned int>::max();

	size_t expectedVerticesDual   = meshFinal.Cell2DsId.size(); // 16 facce triangolate → 16 vertici nel duale
    size_t expectedEdgesDual      = meshFinal.Cell1DsId.size(); // Calcolato dinamicamente, può essere 24 per subdivisionLevel=2
   
    // numero di vertici (id)
	// numero di lati(id)
   
    EXPECT_EQ(meshDual.Cell0DsId.size(), expectedVerticesDual);
    EXPECT_EQ(meshDual.Cell1DsId.size(), expectedEdgesDual);
	
	// ogni vertice del poliedro originale genera una faccia nel duale
	EXPECT_GE(meshDual.Cell2DsId.size(), 4);  // Minimo 4 se è un tetraedro chiuso
	
	// verifico se ogni faccia duale ha almento 3 vertici
	for (const auto& dualFace : meshDual.Cell2DsVertices) {
    	EXPECT_GE(dualFace.size(), 3);
    }
	
	// controllo che uno spigolo non connetta un vertice a se stesso
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		int a = meshDual.Cell1DsExtrema(i, 0);
		int b = meshDual.Cell1DsExtrema(i, 1);
		EXPECT_NE(a, b); // Mi verifica che a e b siano diversi 
    }
	
	// verifico che il calcolo del baricentro sia corretto
	// Per ogni faccia della mesh triangolata
    for (size_t faceId = 0; faceId < meshFinal.Cell2DsId.size(); ++faceId) {
        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
		
        // accedo alla lista dei vertici che compongono la faccia
        const vector<unsigned int>& faceVertices = meshFinal.Cell2DsVertices[faceId];

        // sommo le coordinate dei vertici della faccia
        for (unsigned int v_id : faceVertices) {
			sumX += meshFinal.Cell0DsCoordinates(0, v_id); // coordinata x del vertice v_idMore actions
            sumY += meshFinal.Cell0DsCoordinates(1, v_id); // coordinata y del vertice v_id
            sumZ += meshFinal.Cell0DsCoordinates(2, v_id); // coordinata z del vertice v_id
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

        EXPECT_NEAR(baryX, X, eps);
        EXPECT_NEAR(baryY, Y, eps);
        EXPECT_NEAR(baryZ, Z, eps);
    }

    // verifico che i lati siano costruiti nel modo corretto ---> devono connettere i baricentri
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
	
	// verifico che i lati e le facce siano ordinate correttamente 
	
	for (size_t f = 0; f < meshDual.Cell2DsId.size(); ++f) {
		const auto& edges = meshDual.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshDual.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
       
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (currentEdge == meshDual.Cell1DsId[i]) {
                    currentEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
            
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (nextEdge == meshDual.Cell1DsId[i]) {
                    nextEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext);

            // Recupera il vertice della faccia
            unsigned int faceVertex = vertices[e]; // ID master del vertice della faccia

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


//*********** test triangolazione due



TEST(TestPolyedra, TestOrderedEdges2)
{ 
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generateTetrahedron(mesh);
	
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh meshFinal;
	PolyhedralMesh meshTriangulated2;
	
	int q = 3;
	int b = 2;
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dimension = ComputePolyhedronVEF(q, b, 0);
	vector<int> dimensionDuplicated = CalculateDuplicated(q, b, 0, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, 0,  dimensionDuplicated);
	RemoveDuplicatedVertices(meshTriangulated);
	RemoveDuplicatedEdges(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	
	vector<int> dimension2 = CalculateDimension2(b, q);
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
	triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
	
	// per ogni faccia i lati sono ordinati in modo che la fine dell'arco e coincida con l'inizio dell'arco successivo (e+1)%E
	// il vertice e della faccia deve corrispondere all'origine dell'arco e

    // ciclo su tutte le facce della mesh triangolata 
    for (size_t f = 0; f < meshFinal.Cell2DsId.size(); ++f) {
		const auto& edges = meshFinal.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshFinal.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
        
		
		ASSERT_EQ(vertices.size(), E) << "Numero di vertici e di lati non corrispondono per faccia " << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshFinal.Cell1DsId.size(); i++) {
                if (currentEdge == meshFinal.Cell1DsId[i]) {
                    currentEdgeOrigin = meshFinal.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshFinal.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshFinal.Cell1DsId.size(); i++) {
                if (nextEdge == meshFinal.Cell1DsId[i]) {
                    nextEdgeOrigin = meshFinal.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshFinal.Cell1DsExtrema(i, 1);
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

  
  
  
  
//******************* test cammini minimi (shortest path) ***************************************
//******adapt*******

TEST(Polyhedral_Test, ShortestPath){
	
	PolyhedralMesh mesh;
	PolyhedralMesh meshTriangolazione; //fissata expected
	PolyhedralMesh meshFinal; //in teoria questa è quella "pulita" da scorrere e vedere se abbiamo trovato i shortest paths
	
	//chiamo generazione tetrahedron
	
	generateTetrahedron(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	int startVertexId = 0;
	int endVertexId = 7;	
	
	vector<int> dim = VEF_Calc(q, b, c);
	vector<int> dim_duplicates = duplicates(q, b, c, dimension);
	
//************************

	triangulateAndStore(mesh, meshTriangolazione, b, c, dim_duplicates);
	RemoveDuplicatedEdges(meshTriangolazione);
	RemoveDuplicatedVertices(meshTriangolazione);
	NewMesh(meshTriangolazione, meshFinal, dim);
	
	ShortestPathResult expected(2, 1.63299);
	MatrixXi adjMatrix = calculateAdjacencyMatrix(meshFinal);
	ShortestPathResult result = findShortestPathDijkstra(meshFinal, adjMatrix, startVertexId, endVertexId);
	
	EXPECT_NEAR(result.numEdges, expected.numEdges, 1e-5) 
		<< "error, not the shortest path"
		<< " Expected: " << expected.numEdges << ", Final: " << result.numEdges;
		
	EXPECT_NEAR(result.totalLength, expected.totalLength, 1e-5) 
        << "error, not the shortest path"
        << " Expected: " << expected.totalLength << ", Final: " << result.totalLength;
	 
	}
 
}





//***************************funzioni







/* ************** questo dobbiamo metterlo da un'altra parte ma le teniamo qui per adesso mentre non //fatte push 





namespace PolyhedralLibrary
{
	void inv_val(int& a, int& b) {
		int temp = a;
		a = b;
		b = temp;
	}
	
	vector<int> VEF_Calc(const int q, const int b, const int c)
	{
		vector<int> output(3); // inizializza un vettore con 3 valori, tutti -1 
	
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
	
		output[0] = V;  
		output[1] = E;  
		output[2] = F; 
	
		return output;  // V,E,F
	}






//****************** duplicates

vector<int> duplicates(const int q, const int b, const int c, const vector<int>& dim)
   {
		vector<int> output(3);
		int sub_level = 0;
		sub_level = b + c;
		int V = dim[0];
		int E = dim[1];
		
		if (q == 3) {
			V += 2*4 + 6*(sub_level - 1);
			E += 6*sub_level;
		}
		else if (q == 4) {
			V += 3*6 + 12*(sub_level - 1);
			E += 12*sub_level;
		}
		else {
			V += 4*12 + 30*(sub_level - 1);
			E += 30* sub_level;
		}
		output[0] = V;  
		output[1] = E;  
		output[2] = dim[2]; //il numero delle facce non cambia
		
		return output;
	}

*/