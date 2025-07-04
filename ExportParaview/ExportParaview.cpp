#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>

using namespace std;
namespace PolyhedralLibrary{
	
void ExportParaview(const PolyhedralMesh& meshOutput){
	Gedim::UCDUtilities utilities;
	
	Eigen::VectorXi vert_id(meshOutput.Cell0DsId.size()); //da size_t a Eigen::VectorXi
	for (size_t i = 0; i < meshOutput.Cell0DsId.size(); ++i) //scorriamo 
		vert_id[i] = static_cast<int>(meshOutput.Cell0DsId[i]);
		
	Eigen::VectorXi edge_id(meshOutput.Cell1DsId.size());
	for (size_t i = 0; i < meshOutput.Cell1DsId.size(); ++i)
		edge_id[i] = static_cast<int>(meshOutput.Cell1DsId[i]);
	
   bool markers_Cell0Ds = !meshOutput.Cell0DsMarker.empty();

   bool markers_Cells1Ds = !meshOutput.Cell1DsMarker.empty();

   if (markers_Cell0Ds && markers_Cells1Ds)
   {
	   
	   vector<double> vec_vert;
	   vec_vert.reserve(meshOutput.Cell0DsMarker.size());
	   for (auto m : meshOutput.Cell0DsMarker){
		   vec_vert.push_back(static_cast<double>(m));
		  }
		  Gedim::UCDProperty<double> vert_property;
		  vert_property.Label = "PathVertex";
		  vert_property.UnitLabel = "-";
		  vert_property.NumComponents = 1;
		  
		  vert_property.Data = vec_vert.data();
		  vert_property.Size = vec_vert.size();
		  
		  vector<double> vec_edge;
		  vec_edge.reserve(meshOutput.Cell1DsMarker.size());
		  for (auto m : meshOutput.Cell1DsMarker){
			  vec_edge.push_back(static_cast<double>(m));
			  }
		Gedim::UCDProperty<double> edge_property;
		edge_property.Label = "PathEdge";
		edge_property.UnitLabel = "-";
		edge_property.NumComponents = 1;
		edge_property.Data = vec_edge.data();
		edge_property.Size = vec_edge.size();
		
		utilities.ExportPoints("./Cell0Ds.inp",meshOutput.Cell0DsCoordinates,{ vert_property }, vert_id);
		
		utilities.ExportSegments("./Cell1Ds.inp",meshOutput.Cell0DsCoordinates,meshOutput.Cell1DsExtrema.transpose(),{ vert_property },{ edge_property },edge_id);
	}else{// Export senza proprietà scalari
			utilities.ExportPoints("./Cell0Ds.inp",meshOutput.Cell0DsCoordinates,{},Eigen::VectorXi()); //vuoto se non ha delle propietà
			utilities.ExportSegments("./Cell1Ds.inp",meshOutput.Cell0DsCoordinates,meshOutput.Cell1DsExtrema.transpose(),{},{},Eigen::VectorXi());
	}
	
}

//ci serve ancora?

void MeshTriOutput(const PolyhedralMesh& mesh) {
	//a colonne
		cout << "Cell0DsId: "; 
		for (auto id : mesh.Cell0DsId) cout << id << " ";
		cout << "\nCell0DsCoordinates: " << endl;
		for (int j = 0; j < mesh.Cell0DsCoordinates.cols(); ++j) {
			cout << "column " << j << ": ";
			for (int i = 0; i < mesh.Cell0DsCoordinates.rows(); ++i) {
				cout << mesh.Cell0DsCoordinates(i, j) << " ";
			}
			cout << endl;
		}
		cout << "Cell0DsFlag: " << endl;
		for (const auto& row : mesh.Cell0DsFlag) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		
		//lati
		
		//a righe
		cout << "Cell1DsId: "; 
		for (auto id : mesh.Cell1DsId) cout << id << " ";
		cout << "\nCell1DsExtrema: " << endl;
		for (int i = 0; i < mesh.Cell1DsExtrema.rows(); ++i) {
			for (int j = 0; j < mesh.Cell1DsExtrema.cols(); ++j) {
				cout << mesh.Cell1DsExtrema(i,j) << " ";
			}
			cout << endl;
		}
		cout << "Cell1DsFlag:" << endl;
		for (const auto& row : mesh.Cell1DsFlag) {
			cout << row << " ";
			cout << endl;
		}
		cout << "Cell1DsOriginalFlag:" << endl;
		for (const auto& row : mesh.Cell1DsOriginalFlag) {
			cout << row << " ";
			cout << endl;
		}
		
		//facce
		cout << "Cell2DsId: "; 
		for (auto id : mesh.Cell2DsId) cout << id << " ";
		
		cout << "\nCell2DsVertices: " << endl;
		for (const auto& row : mesh.Cell2DsVertices) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		cout << "Cell2DsEdges: " << endl;
		for (const auto& row : mesh.Cell2DsEdges) {
			for (auto v : row) cout << v << " ";
			cout << endl;
		}
		
		cout <<"Cell3DsId: " << mesh.Cell3DsId << endl; 
		cout <<"nº vertici: " << mesh.NumCells0Ds << endl;
		cout <<"nº lati: " << mesh.NumCells1Ds << endl;
		cout <<"nº facce: " << mesh.NumCells2Ds << endl;
		
		cout <<"Cell3DsV Vertici: "; 
		for (auto id : mesh.Cell3DsVertices) cout << id << " ";
		cout << endl;
		cout <<"Cell3Ds Lati: "; 
		for (auto id : mesh.Cell3DsEdges) cout << id << " ";
		cout << endl;
		cout <<"Cell3Ds Faccie: "; 
		for (auto id : mesh.Cell3DsFaces) cout << id << " ";
		cout << endl;

		cout <<"\n"<< endl;
	}
}