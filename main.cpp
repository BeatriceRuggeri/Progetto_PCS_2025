int main()
{
	bool dual_flag;
    int P;
	int Q;
	int B;
	int C;
	int id_vertice1;
	int id_vertice2;

    cout << "Inserisci un valore per P: ";
    cin >> P;
	
	if (P!=3) {
		cerr << "Valore di P non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per P: ";
		cin >> P;
	}
	else {
		cout << "Hai inserito: " << P << endl;
	}
	
	cout << "Inserisci un valore per Q: ";
    cin >> Q;
	
	if (Q!=3 & Q!= 4 & Q!=5) {
		cerr << "Valore di Q non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per Q: ";
		cin >> Q;
	}
	else {
		cout << "Hai inserito: " << Q << endl;
	}
	
	cout << "Inserisci un valore per B: ";
    cin >> B;
	if (B < 1) {
		cerr << "Valore di B non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per B: ";
		cin >> B;
	}
	cout << "Inserisci un valore per C: ";
    cin >> C;
	if (C < 1) {
		cerr << "Valore di C non valido, si prega di riprovare" << endl;
		cout << "Inserisci un valore per C: ";
		cin >> C;
	}
	
	vector<int> Quadrupla = {P, Q, B, C} 
	cout << Quadrupla << endl;
	

	PolyhedralMesh mesh;
	PolyhedralMehs meshOutput;

	
    
 
	if (P == 3 && B != C){
		dual_flag = false;
		if (Q == 3){
			gen_tetraedro(mesh);
		}
		if (Q == 4){
			gen_ottaedro(mesh);
		}
		if (Q == 5){
			gen_icosaedro(mesh);
		}

		int new_id_v = mesh.M0D.size() + 1;
    	int new_id_s = 1;
    	int new_id_f = 1;

		triangolazione1(B, new_id_v, new_id_s, new_id_f, mesh.M0D, mesh.M1D, mesh.M2D, meshOutput.M0D, meshOutput.M1D, meshOutput.M2D);

		vector<vector<double>> tabella_confinamenti; //tab_conf
		vector<vector<double>> vertici_duale; //v_duale
		vector<vector<double>> spigoli_duale; //s_duale

		if (dual_flag == true){
			costruttore_tabella_confinamenti_3(meshOutput.M2D,tabella_confinamenti);
			baricentro_3(meshOutput.M0D, meshOutput.M2D, vertici_duale);
			costruttore_duale_3(tabella_confinamenti, meshOutput.M0D, meshOutput.M2D, vertici_duale, spigoli_duale);
		}
	} else if (P == 4 && B != C){
		costruttore_tabella_confinamenti_4(meshOutput.M2D,tabella_confinamenti);
		baricentro_4(meshOutput.M0D, meshOutput.M2D, vertici_duale);
		costruttore_duale_4(tabella_confinamenti, meshOutput.M0D, meshOutput.M2D, vertici_duale, spigoli_duale);
	}
	

	


	GedimUCDUtilities utilities;
    utilities.ExportPoints("./Cell0Ds.inp",
                            mesh.M0D);
    utilities.ExportSegments("./Cell1Ds.inp",
                              mesh.M0D,
                              mesh.M1D);
                              
	
}
	
	
	/*
	int Matrix PQ[5][2] = {{3, 3}, {3, 4}, {4, 3}, {5, 3}, {3, 5}};
	
	for (i= 0; i<=5; i++)
	{	
		if (PQ[i]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE1(P, Q, b, c)
		}	
		if (PQ[i][1]==3)
		{
			int P = PQ[i][0];
			int Q = PQ[i][1];
			// CLASSE2 (P, Q, b, c)
		}
	}
	*/
	//trial
	
	PolyhedralLibraryPolyhedralMesh mesh;
	/*
	bool success_edges_cube = ExportCell1Ds("C_Cell1Ds.txt",mesh.cube_edges,12);
	
	if(!success_edges_cube){
		cerr<<"Export cube's edges failed."<<endl;
		return 1;
	}
	
	bool success_tetrahedron = ExportCell0Ds("T_Cell0Ds.txt",mesh.Vert_tetrahedron,4);
	
	if(!success_tetrahedron){
		cerr<<"Export tetrahedron failed."<<endl;
		return 1;
	}
	
	
	bool success_cube = ExportCell0Ds("C_Cell0Ds.txt",mesh.Vert_cube,8);
	
	if(!success_cube){
		cerr<<"Export cube failed."<<endl;
		return 1;
	}
	
	bool success_octahedron = ExportCell0Ds("O_Cell0Ds.txt",mesh.Vert_octahedron,6);
	
	if(!success_octahedron){
		cerr<<"Export octahedron failed."<<endl;
		return 1;
	}
	
	bool success_dodecahedron = ExportCell0Ds("D_Cell0Ds.txt",mesh.Vert_dodecahedron,20);
	
	if(!success_dodecahedron){
		cerr<<"Export dodecahedron failed."<<endl;
		return 1;
	}
	
	bool success_icosahedron = ExportCell0Ds("I_Cell0Ds.txt",mesh.Vert_icosahedron,12);
	
	if(!success_icosahedron){
		cerr<<"Export icosahedron failed."<<endl;
		return 1;
	}
	*/			