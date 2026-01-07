/**************************
 * Copyright (©) 2025 FAIR AERO. All rights reserved.
 * This material is proprietary and confidential to FAIR AERO.
 * Its use and distribution are subject to the terms of a Non-Disclosure Agreement (NDA)
 * and are authorized only for educational purposes by designated participants.
 * This script forms part of the GENEPI™ software, a trademark of FAIR AERO.
 * For further information, please contact: info@fair.aero
 **************************/

#include "structuralSolver.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

void StructuralSolver::loadMesh(string meshfile) {
	///Read mesh file from path.
	string line; 
	ifstream mfile(meshfile);
	if (mfile.is_open()) {
		getline(mfile, line);
		vector<string> cut = Utils::split(line,',');
		cout << "reading " << cut[1] << " points... ";
		int np = stoi(cut[1]);
		nPoints = np;
		getline(mfile, line);
   	
		//Read points
		for(int i=0; i< np; i++) {
			vector<string> linecut;
			Utils::tokenize(line,linecut,",");
			double x = stod(linecut[0]);
			double y = stod(linecut[1]);
			double z = stod(linecut[2]);
			Array vertex(x,y,z);
			Point* pt = new Point(vertex);
			pt->Xnew = vertex;
			pt->index = (int)i;		
			pt->v = Array(0.,0.,0.); 
			pt->vnew = Array(0.,0.,0.);
			pt->fex = Array(0.,0.,0.);
			pt->res = 0.0;
			pt->phi = Array(0.,0.,0.);
			pt->err = 0.;
			points.push_back(pt);
			getline(mfile, line);
		}

		//Read facets
		cut = Utils::split(line,',');
		nFacets = stoi(cut[1]);
		getline(mfile, line);
		cout << nFacets << " facets ... "; 
		facets.clear();	
		vector<Facet*> sliver_facets;

		double minAR = 1.;
		double avgAR = 0.;
		double avgAR2 = 0.; 
		double maxAR = 1.;
		int nSliver = 0; 
		int nBadQuality = 0; 

		for(int i=0; i< nFacets ; i++) {
			vector<string> cut4;
			Utils::tokenize(line,cut4,",");
			size_t id1 = stoi(cut4[0]);
			size_t id2 = stoi(cut4[1]);
			size_t id3 = stoi(cut4[2]);

			if(id1 >= points.size() || id2 >= points.size() || id3 >= points.size()) {
				cout << "Error with mesh file for facet " << id1 << " " << id2 << " " << id3 << endl; 
				cout << "Number of points " << points.size() << endl; 
				exit(1); 
			}

			Point* p1 = points[id1]; 
			Point* p2 = points[id2]; 
			Point* p3 = points[id3]; 

			Facet * facet = new Facet(p1,p2,p3);
			facet->sigma_1 = 0.; 
			facet->sigma_2 = 0.; 
			facet->sigma_vm = 0.;
			facet->setYoung(Young,poisson);
			
			double AR = facet->getAR(); 
			if(AR > maxAR) maxAR = AR; 	
			if(AR < minAR) minAR = AR;
			avgAR += AR; 
			if(facet->isFuckedUp()) {
				sliver_facets.push_back(facet);
				nSliver ++; 
			}
			else if(AR > ARcrit) {
				sliver_facets.push_back(facet);
				nBadQuality++; 
			}	
			else {
				this->facets.push_back(facet); 
				avgAR2 += AR; 
			}

			getline(mfile, line);
		}
		avgAR /= (double)nFacets; 
		avgAR2 /= (double)facets.size(); 	

		if(sliver_facets.size() > 0) {
			string sliverfile = path+"sliver_triangles.stl"; 
			Utils::writeToStl(sliverfile,sliver_facets); 
		}

        // Write mesh quality report
		cout << endl; 
		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl; 
		cout << "%%%%%%%%%%% MESH QUALITY REPORT %%%%%%%%%%%%%%%%%%%%%%%%%" << endl; 
		cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl; 
		cout << endl; 
		cout << "%% Average facet AR (before cleaning) = " << avgAR << endl;
	    cout << "%% Average facet AR (after cleaning) = " << avgAR2 << endl; 
		cout << "%% Max AR = " << maxAR << " with a cutoff at ARcrit = " << ARcrit << endl; 
		cout << "%% Min AR = " << minAR << endl; 
		cout << "%% Deleted " << nSliver << " sliver facets." << endl; 
		cout << "%% Deleted " << nBadQuality << " facets with poor quality." << endl;
	    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl; 
		cout << endl; 	

		mfile.close();

	}else {
		cout << "Could not open mesh file: " << meshfile << endl; 
		exit(1);
	}

    //Compute area & masses
	computeMassesAndAreas();

	cout << "Mesh loaded." << endl; 
}

void StructuralSolver::rotate(double angle_deg,Array direction,Array center) {
	Eigen::Matrix<double,3,3> rotmat(3,3); 
	
	direction.normalize();
	double angle_rad = EIGEN_PI * angle_deg / 180.0;
	
	double c = cos(angle_rad);
	double s = sin(angle_rad);
	
	double ux = direction[0];
	double uy = direction[1];
	double uz = direction[2];

	rotmat.coeffRef(0, 0) = ux * ux*(1 - c) + c;
	rotmat.coeffRef(0, 1) = ux * uy*(1 - c) - (uz*s);
	rotmat.coeffRef(0, 2) = ux * uz*(1 - c) + (uy*s);

	rotmat.coeffRef(1, 0) = ux * uy*(1 - c) + (uz*s);
	rotmat.coeffRef(1, 1) = uy * uy*(1 - c) + c;
	rotmat.coeffRef(1, 2) = uy * uz*(1 - c) - (ux*s);

	rotmat.coeffRef(2, 0) = ux * uz*(1 - c) - (uy*s);
	rotmat.coeffRef(2, 1) = uy * uz*(1 - c) + (ux*s);
	rotmat.coeffRef(2, 2) = uz * uz*(1 - c) + c;

	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i];
		//translation pour que le centre soit à l'origine du repère lors de la rot
		Array rotated = rotmat * (p->X - center) + center;
		p->X = rotated; 
		p->Xnew = rotated;
	}	
}

void StructuralSolver::translate(Array translation) {
	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i]; 
		p->X += translation; 
		p->Xnew = p->X; //ADDED from Tom
	}
}

void StructuralSolver::initializeDelingette() {
	//Build an array of edges, based on the triangulation.
	//each edge must be uniquely defined, this is why we start
	//with a map.
	cout << "Build edges for triangle based simulation ... "; 
	map<pair<int,int>,vector<Facet*>> edge_map; 	
	for(size_t i=0; i<facets.size(); i++) {
		Facet* facet = facets[i]; 

		//3 edges of the triangle
		pair<int,int> edge1,edge2,edge3;
		//Edge 1 
		int ind0 = facet->points[0]->index; 
		int ind1 = facet->points[1]->index;
		if(ind0 < ind1) {
			edge1.first = ind0; edge1.second = ind1;
		} 
		else {
			edge1.first = ind1; edge1.second = ind0; 
		}

		//Edge 2 
		int ind2 = facet->points[2]->index;
		if(ind1 < ind2) {
			edge2.first = ind1; edge2.second = ind2;
		}
		else {
			edge2.first = ind2; edge2.second = ind1;
		}
			
		//Edge 3	
		if(ind2 < ind0) {
			edge3.first = ind2; edge3.second = ind0; 
		}
		else {
			edge3.first = ind0; edge3.second = ind2; 
		}

		std::map<pair<int,int>,vector<Facet*>>::iterator it = edge_map.find(edge1);
	       	if(it == edge_map.end()) { 
			vector<Facet*> adjacent_facets; 
			adjacent_facets.push_back(facet);

			edge_map[edge1] = adjacent_facets; 
		}
		else{ //add the triangle 
			it->second.push_back(facet); 
		}
		
		it = edge_map.find(edge2); 
		if(it == edge_map.end()) {
			vector<Facet*> adjacent_facets; 
			adjacent_facets.push_back(facet);

			edge_map[edge2] = adjacent_facets; 
		}
		else { 
			it->second.push_back(facet);
		}

		it = edge_map.find(edge3);
		if(it == edge_map.end()) {
			vector<Facet*> adjacent_facets; 
			adjacent_facets.push_back(facet);

			edge_map[edge3] = adjacent_facets; 	
		}
		else {
			it->second.push_back(facet); 
		}

	}
	map<pair<int,int>, int> key2index;
	int i = 0;
	for (pair<pair<int,int>,vector<Facet*>> entry : edge_map) {
		key2index[entry.first] = i++;
	}
	

	//Create the Edge to facet list from the resulting hashmap
	map<pair<int,int>,vector<Facet*>>::iterator it; 
	for(it = edge_map.begin(); it!= edge_map.end(); it++) {
		pair<int,int> edge_indexes = it->first;
	       	vector<Facet*> adjacent_triangles = it->second; 

		//Create local edge and compute its geometry.	
		Edge* edge = new Edge(points[edge_indexes.first],points[edge_indexes.second]);		
	
		edge->cvector.resize(0); 

		if(edge->points.first->mass < 1.e-19 || edge->points.second->mass < 1.e-19) {
			cout << "Error one end of an edge has zero mass " << endl; 
			exit(1); 
		}
		edge->computeRestLength(); 
		edges.push_back(edge);
		
		edge->computeLength();
		for(size_t i=0; i< adjacent_triangles.size(); i++) {	
			Facet* facet = adjacent_triangles[i]; 
			facet->edges.push_back(edge);	
		}
	}
	key2index.clear();

	//Compute the stiffnesses 
	for(size_t i=0; i<facets.size(); i++) {
		Facet* facet = facets[i]; 

		assert(facet->edges.size()==3); 

		facet->edge1 = facet->edges[0]; 
		facet->edge2 = facet->edges[1]; 
		facet->edge3 = facet->edges[2]; 
	
		Edge* edge1 = facet->edge1;
		Edge* edge2 = facet->edge2;
		Edge* edge3 = facet->edge3;

		facet->edges.clear();

		//adjust the numbering of the points to stay within Delingette's 
		//framework.
		facet->points.clear();
		
		//Find Q1: intersection between Edge2 and Edge3  
		if(edge2->points.first->index == edge3->points.first->index) facet->points.push_back(edge2->points.first);
		if(edge2->points.first->index == edge3->points.second->index) facet->points.push_back(edge2->points.first);
		if(edge2->points.second->index == edge3->points.first->index) facet->points.push_back(edge2->points.second);
		if(edge2->points.second->index == edge3->points.second->index) facet->points.push_back(edge2->points.second);
		
		//Find Q2: intersection between Edge1 and Edge3
		if(edge1->points.first->index == edge3->points.first->index) facet->points.push_back(edge1->points.first);
		if(edge1->points.first->index == edge3->points.second->index) facet->points.push_back(edge1->points.first);
		if(edge1->points.second->index == edge3->points.first->index) facet->points.push_back(edge1->points.second);
		if(edge1->points.second->index == edge3->points.second->index) facet->points.push_back(edge1->points.second);
	
		//Find Q3: intersection between Edge1 and Edge2
		if(edge1->points.first->index == edge2->points.first->index) facet->points.push_back(edge1->points.first);
		if(edge1->points.first->index == edge2->points.second->index) facet->points.push_back(edge1->points.first);
		if(edge1->points.second->index == edge2->points.first->index) facet->points.push_back(edge1->points.second);
		if(edge1->points.second->index == edge2->points.second->index) facet->points.push_back(edge1->points.second);
	
		assert(facet->points.size()==3);
		assert(facet->points[0]->index!=facet->points[1]->index);
		assert(facet->points[1]->index!=facet->points[2]->index);
		assert(facet->points[2]->index!=facet->points[0]->index);
	}
	

	//Direction of the edges
	for(size_t i=0; i<facets.size(); i++) {
			Facet* facet = facets[i]; 
			facet->computeInitialLength();//check if only once

			//Creating an array to have the direction of each edge 
	 		facet->directions = Array(0.,0.,0.);
		
			//Checking the direction of edgei->dX
			Matrix3 sumX;
       			for(int i=0; i<3; i++){
			
				sumX.coeffRef(i,0)= facet->dX1[i] + facet->edge1->dX[i];
				sumX.coeffRef(i,1)= facet->dX2[i] + facet->edge2->dX[i];
				sumX.coeffRef(i,2)= facet->dX3[i] + facet->edge3->dX[i];
				
			}
			
			//Loop on edges of the facet
			for(int i=0; i<3; i++){
				//First case : sumX1 coefficients are nulls -> dX1 is opposed to edge1->dX
				if ((abs(sumX.coeffRef(0,i)) + abs(sumX.coeffRef(1,i))  + abs(sumX.coeffRef(2,i)))<=1e-6){ 
		
					facet->directions[i] = -1;
					
				}
				//Second case : dX1 is edge1->dX
				else {
					facet->directions[i] = 1;
					
				}
			}
		}

	for(size_t i=0; i<facets.size(); i++) {
		Facet* facet = facets[i];
		
		facet->initializeTRQS();

		//identify points that are connected to facets.
		facet->points[0]->inFacet = true; 
		facet->points[1]->inFacet = true; 
		facet->points[2]->inFacet = true; 
	}

	//Reset phi vectors just in case
	for(size_t i=0; i<points.size(); i++) points[i]->phivector.clear(); 
	
	//Compute the damping coefficients and initPhiVectors
	double avgk = 0.; 
	for(size_t i=0; i<edges.size(); i++) {
		Edge* edge = edges[i]; 
		Point* p1 = edge->points.first;
		Point* p2 = edge->points.second;
	
		//Computing k, order of magnitude
		if (argyris) edge->k = Young;
				
		//Adding sigma damping
		if(sigma<0.) {
			edge->damp_coef = sigma * sqrt(edge->k * (p1->mass + p2->mass));// / edge->getRestLength();  	
		} else { edge->damp_coef = 0.; }
		
		//Adding ksi damping
		if(ksi<0.) {
			edge->damp_coef2 = ksi * sqrt(edge->k * (p1->mass + p2->mass)); //sqrt commented  	
		} else { edge->damp_coef2 = 0.;}

		if (delingette){
			edge->initPhiVector(); 
		}

		//Length scale.
		avgk += edge->k; 
	}
	avgk /= (double)edges.size(); 

	//setup beta damping. 
	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i]; 
		if(useScaledBetaDamping) p->beta = beta_damp * sqrt(avgk * p->mass);
		else p->beta = beta_damp; 
	}

	cout << "Done." << endl; 
}

void StructuralSolver::setOrthotropicConstants(double young_ext_1, double poisson_ext_12, double young_ext_2, double poisson_ext_21,double G_ext_12) {
	double inv_poisson = 1./(1.-poisson_ext_12*poisson_ext_21);
	double coeff1 = young_ext_1 * inv_poisson;
	double coeff2 = young_ext_2 * inv_poisson;
	double coeff12 = poisson_ext_21 * young_ext_1*inv_poisson;
	double coeff21 = poisson_ext_12* young_ext_2* inv_poisson;
	assert(abs(coeff12-coeff21) <=1e-3);
	
	// Filling the matrix with Argyris-Tenek paper
	// Strain-stress matrix along cloth filled with material data
	
	MatData.coeffRef(0,0) = coeff1;
	MatData.coeffRef(0,1) = coeff12;
	MatData.coeffRef(0,2) = 0.;
	MatData.coeffRef(1,0) = coeff21;
	MatData.coeffRef(1,1) = coeff2;
	MatData.coeffRef(1,2) = 0.;
	MatData.coeffRef(2,0) = 0.;
	MatData.coeffRef(2,1) = 0.;
	MatData.coeffRef(2,2) = G_ext_12;//We get Gammaxy in function of tauxy
	
	cout <<"MatData in setOrthotropicConstants  : " <<MatData<< endl;


	for(size_t i=0; i<facets.size(); i++) {
		Facet* facet = facets[i];
		facet->MatDataFacet = MatData;
		facet->Young_weft = young_ext_1;
		facet->Young_warp = young_ext_2;	
	}

}

void StructuralSolver::initializeArgyris(){
	size_t j = 0;

	//Reset PhiVector in case
	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i]; 
		p->phivector.clear(); 
		p->jacvector.clear(); 
		p->jacNonDiagTermArgy.clear();
	}
	for(size_t i=0; i<edges.size(); i++) {
		Edge* edge = edges[i]; 
		edge->jacvector.clear(); 
	}

	
	//Compute Edge Lengthes
	for(j=0;j<edges.size();j++){
		Edge* edge = edges[j];

		edge->computeRestLength();
		edge->computeLength();

		edge->init_jacvector(); 
	}

	//Compute Length Facet and initialize mat axis for flag
	for(j=0; j<facets.size(); j++) {
		Facet* facet = facets[j]; 

		facet->computeInitialLength();
		facet->computeLength(); 
		facet->matAxis = matAxis;
		facet->initPhiVectors();
		facet->init_jacvector_diag();
		facet->init_jacvector_non_diag(); 
		facet->MatDataTemp = facet->MatDataFacet; 
	}
	

	//Beta Damping 
	for(j=0; j<points.size(); j++) {
		Point* p = points[j]; 
		if(useScaledBetaDamping) p->beta = beta_damp * sqrt(Young * p->mass); 	
		else p->beta = beta_damp; 
	}
	
	//Mat Axis
	projectMatAxisOnFacet(); 

	//Compute k_n D_n of Pauletti
	computeNaturalandLocalMatrix();
	
	//Sanity check
	for(j=0; j<facets.size(); j++) {
		Facet* facet = facets[j]; 
		assert(facet->Young_weft > 10.);
		assert(facet->Young_warp > 10.); 
	}

}

void StructuralSolver::computeNaturalandLocalMatrix(){
	size_t i =0;
	#pragma omp parallel private(i)
	{
	#pragma omp for
	for(i=0; i<facets.size(); i++) {
		//Computing Transfo_B0_Bl and Transfo Bc_Bl
		facets[i]->createTransfoMatrix();
		
		//Compute Dn and k_n
		facets[i]->computeNaturalandLocalMatrix();	
	}
	}
}

void StructuralSolver::projectMatAxisOnFacet(){
	size_t i=0;
	#pragma omp parallel private(i)
	{
	#pragma omp for
 
	for(i=0; i<facets.size(); i++) {
		Facet* facet = facets[i];
		facet->createTransfoB0Bl2();

		//Projection of material_axis on the facet
		Array point_material_0 = Array(0.,0.,0.);//Global Cartesien frame
       		Array point_material_1 = point_material_0 + facet->matAxis;//Global Cartesian frame
		Array plane_point = facet->points[0]->Xnew;//Global Cartesian frame
		Array plane_normal = Array(facet->Transfo_B0_Bl(0,2),facet->Transfo_B0_Bl(1,2),facet->Transfo_B0_Bl(2,2));//w_l
			
		Array material_0_proj = point_material_0 - (point_material_0 - plane_point).dot(plane_normal)*plane_normal;//Projection of point_0
		Array material_1_proj = point_material_1 - (point_material_1 - plane_point).dot(plane_normal)*plane_normal;//Projection of point_1
		
		//Material Axis projected but expressed in global frame
		Array matAxisProjected =  material_1_proj - material_0_proj;

		//Material Axis projected in local frame (thus the visualisation is random)
		facet->matAxis_local = facet->Transfo_B0_Bl.inverse()*matAxisProjected;
	}	
	}
}

void StructuralSolver::applyGravity(double gravity_constant) {
	Array total_force(0.,0.,0.); 
	for(size_t i=0; i<points.size(); i++) {
		Point* pi = points[i];
		double mg = 0.; 
		mg = -gravity_constant * pi->mass;
		pi->fex = Array(0.,0.,mg); 
		total_force+=pi->fex;
	}
	cout << " Total gravity force = " << total_force << endl; 
}

/** Writes the entire structure to VTK file
 */
 void StructuralSolver::printVTK(string file) {
    ofstream myfile(file); //Main VTK file. 

    cout << "writing structure to : " << file << endl;

    if (myfile.is_open()){
        myfile << "# vtk DataFile Version 3.0" << endl;
        myfile << "Data generated by GENEPI" << endl;
        myfile << "ASCII" << endl;
        myfile << "DATASET UNSTRUCTURED_GRID" << endl;
        myfile << "POINTS " << this->points.size() << " double" << endl;

        int np = points.size();
        int nf = facets.size();

        for(size_t i=0; i<points.size(); i++) {
        Point* p = points[i]; 
        Array X = p->X;
        double x = X[0];
        double y = X[1];
        double z = X[2];

        //Write x,y,z with high accuracy in case simulation has to be 
        //re-loaded.		
        stringstream ssx;
        ssx << setprecision(10) << x;
        string str_x;
        ssx >> str_x; 

        stringstream ssy;
        ssy << setprecision(10) << y;
        string str_y;
        ssy >> str_y; 

        stringstream ssz;
        ssz << setprecision(10) << z;
        string str_z;
        ssz >> str_z; 

        myfile << str_x << " " << str_y << " " << str_z << endl;
        }
        int ndata = 4 * (int)this->facets.size();
        myfile << "CELLS " << this->facets.size() << " " << ndata << endl;

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 

        int id1 = facet->points[0]->index;
        int id2 = facet->points[1]->index;
        int id3 = facet->points[2]->index;

        myfile << "3 " << to_string(id1) << " " << to_string(id2) << " " << to_string(id3) << endl;;	
        }

        myfile << "CELL_TYPES " << nf << endl;

        for(size_t i=0; i<facets.size(); i++) {
        myfile << "5" << endl;
        }

        myfile << "CELL_DATA " << nf << endl; 

        myfile << "VECTORS stresses double" << endl; 	

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 

        string sep = " ";	
        string topr = to_string(facet->sigma_1) + sep + to_string(facet->sigma_2) + sep + to_string(facet->sigma_vm);

        myfile << topr << endl;
        }

        myfile << "VECTORS matAxis double" << endl; 	

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 

        string sep = " ";	
        string topr = to_string(facet->matAxis[0]) + sep + to_string(facet->matAxis[1]) + sep + to_string(facet->matAxis[2]);

        myfile << topr << endl;
        }

        myfile << "VECTORS normals double" << endl; 	

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 

        string sep = " ";	
        string topr = to_string(facet->normal[0]) + sep + 
        to_string(facet->normal[1]) + sep + to_string(facet->normal[2]);

        myfile << topr << endl;
        }

        myfile << "SCALARS crease_criterion double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl;

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 
        string topr = to_string(facet->crease_crit);

        myfile << topr << endl;
        }

        myfile << "SCALARS svm_max double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl;

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 
        string inrad = to_string(facet->svm_max);

        myfile << inrad << endl;
        }

        myfile << "SCALARS Young_weft double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl;

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 
        string topr = to_string(facet->getYoung()); 

        myfile << topr << endl;
        }

        myfile << "SCALARS Young_warp double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl;

        for(size_t i=0; i<facets.size(); i++) {
        Facet* facet = facets[i]; 
        string Y_warp = to_string(facet->Young_warp); 

        myfile << Y_warp << endl;
        }

        myfile << "POINT_DATA " << np << endl;

        myfile << "SCALARS Mass double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl; 

        for(int i=0; i< np; i++) {
        Point* p = points[i]; 
        string topr =  to_string(1./p->invmass); 

        myfile << topr << endl; 
        }

		myfile << "SCALARS area double 1" << endl; 
		myfile << "LOOKUP_TABLE default" << endl;

		for(int i=0; i< np; i++) {
			Point* p = points[i]; 
			string topr = to_string(p->area);

			myfile << topr << endl; 
		}

        myfile << "SCALARS error double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl; 

        for(int i=0; i< np; i++) {
        Point* p = points[i]; 
        string topr =to_string(p->err);

        myfile << topr << endl; 
        }
        myfile << "SCALARS residuals double 1" << endl; 
        myfile << "LOOKUP_TABLE default" << endl; 

        for(int i=0; i< np; i++) {
        Point* p = points[i]; 
        stringstream ssres;

        ssres << setprecision(10) << p->res;
        string str_res;
        ssres >> str_res; 
        string topr = str_res; 

        myfile << topr << endl; 
        }

        myfile << "VECTORS External_forces double" << endl;

        for(int i=0; i< np; i++) {
        Point* p = points[i];

        double fp_x = p->fex[0];
        double fp_y = p->fex[1];
        double fp_z = p->fex[2];
        string sep = " ";
        string topr = to_string(fp_x) + sep + to_string(fp_y) + sep + to_string(fp_z);

        myfile << topr << endl; 
        }

        myfile << "VECTORS normals double" << endl;

        for(int i=0; i< np; i++) {
        Point* p = points[i];

        double nx = p->normal[0];
        double ny = p->normal[1];
        double nz = p->normal[2];

        string sep = " ";
        string topr = to_string(nx) + sep + to_string(ny) + sep + to_string(nz);

        myfile << topr << endl; 
        }

        myfile.close();
    }else{
        cout << "Unable to open file";
    }
}

void StructuralSolver::advanceToTimeEuler(double current_time, double next_time, double& dt_structure, int& total_struct_iters) {
    //Advance to time perfoming the necessary number of solver iterations.
	double time = current_time;
	double current_dt = dt_structure;  
	bool success = false;		

	cout << "Starting advanceToTime with current time : " << current_time << " s going to next time = " << next_time << " s." << endl;  

	while(time < next_time) {
		
		double future_time = time + dt_structure; 
		if((future_time > next_time) && use_variable_timestep) {
			double dt_structure_ = next_time - time; 
			if(dt_structure_ < 1.e-12) { 
				time = next_time; 
			}
			else {
				dt_structure = dt_structure_; 
			}
		}
		if(time==next_time) break; //This is to avoid tiny timestep if there is a zero error while finishing the iteration 

		current_dt = dt_structure;
		success = advanceEulerImplicitBaraffWitkin(use_variable_timestep,dt_structure);

		if(success) {
			time += current_dt; 
			total_struct_iters++;
			cout << "accepted step " << total_struct_iters << " at t = " << time << " with err = " 
				<< err << " and res = " << res << " with dt = " << current_dt << " next dt = " << dt_structure << endl;
		
			if(current_dt < 1.e-15) {
				cout << "Simulation killed after dt went below 1e-15. That sucks ... " << endl; 
				exit(1); 
			}
		}	
		else {
			cout << " REJECTED step, using dt = " << dt_structure << endl; 
			current_dt = dt_structure;
			success = advanceEulerImplicitBaraffWitkin(use_variable_timestep,dt_structure);
			while(!success) {
				cout << " REJECTED step, using dt = " << dt_structure << endl; 
				success = advanceEulerImplicitBaraffWitkin(use_variable_timestep,dt_structure);
				current_dt = dt_structure; 	
			}
			time += current_dt;
			total_struct_iters++;
			cout << "accepted step " << total_struct_iters << " at t = " << time << " with err = " 
				<< err << " and res = " << res << " with dt = " << current_dt << " next dt = " << dt_structure << endl;
			
			if(current_dt < 1.e-15) {
				cout << "Simulation killed after dt went below 1e-15. That sucks ... " << endl; 
				exit(1); 
			}		
		}
	}
	double time_error = time - next_time ; 
	cout << "time error = " << time_error <<  endl; 
}

bool StructuralSolver::advanceEulerImplicitBaraffWitkin(bool use_variable_timestep, double& dt_) {
	bool printThings = false; 

	int n = points.size(); 
	bool accepted = true; 
	int max_iters = max_solver_it;

	Eigen::setNbThreads(num_threads);
	omp_set_num_threads(num_threads);
	
	spm jX(3*n,3*n); //, Eigen::RowMajor	
	spm jV(3*n,3*n); //, Eigen::RowMajor	

	bool iterative = true; 

	Eigen::DGMRES<spm, Eigen::IncompleteLUT<double> > solver_it;

	if(iterative) {
		solver_it.preconditioner().setDroptol(drop_tol);
		solver_it.preconditioner().setFillfactor(fill_factor);
		if(use_variable_timestep) solver_it.setMaxIterations(max_iters);
	}

	double time = omp_get_wtime();
					
	//compute sparse identity matrix
	spm identity(3*n,3*n);
	vector<T> jAl; 
	for(int i=0; i<3*n; i++) {
		jAl.push_back(T(i,i,1.));		
	}
	identity.setFromTriplets(jAl.begin(),jAl.end());	
	identity.makeCompressed();
	
	//compute sparse inverse mass matrix
	spm invmass(3*n,3*n); //, Eigen::RowMajor	
	jAl.clear(); 
	for(int i=0; i<n; i++) {
		double im = points[i]->invmass; 	
		jAl.push_back(T(3*i,3*i,im));
		jAl.push_back(T(3*i+1,3*i+1,im));
		jAl.push_back(T(3*i+2,3*i+2,im));		
	}
	invmass.setFromTriplets(jAl.begin(),jAl.end());	
	invmass.makeCompressed();
	
	//Solve with deltaV baraff & Witkin style	

	//Initialisation for the first time step. 
	//After that we use the old dv value as an initial guess.
	//Unless the time step has been rejected.  
	if(bw_init) {
		dv.resize(3*n);
		for(int i=0;i<3*n;i++) {
			dv[i] = 0.0; 
		}
		bw_init = false; 

		cout << "dv vector reset" << endl; 
	}
	
	time = omp_get_wtime() - time;
	
	if(printThings)	cout << "initialization step : " << time << endl; 
	
	//Compute jacobian and internal forces. 
	updateInternalForces();

	time = omp_get_wtime(); 
	computeJacobian(jX,jV);		
	
	time = omp_get_wtime() - time;	
	
	if(printThings){
		cout << "Jacobian step : " << time << ", jxnorm : " << jX.norm() << ", jvnorm : " << jV.norm() << endl; 
	}

	Eigen::VectorXd v_zero(3*n); 
	Eigen::VectorXd phi_zero(3*n); 
	for(int i=0; i<n; i++) {
		Point* pi = points[i]; 
		Array forces = pi->phi + pi->fex; 
		phi_zero.coeffRef(3*i) =   forces[0]; 
		phi_zero.coeffRef(3*i+1) = forces[1]; 
		phi_zero.coeffRef(3*i+2) = forces[2]; 
		v_zero.coeffRef(3*i) =   pi->v[0]; 
		v_zero.coeffRef(3*i+1) = pi->v[1]; 
		v_zero.coeffRef(3*i+2) = pi->v[2]; 

	}

	//compute LHS, for weird compile reasons we have to
	//add identity in a second line
	spm LHS(3*n,3*n); //, Eigen::RowMajor	
	LHS = - dt_*invmass*(jV + dt_*jX);//+ dt*jV;//(identity  - dt*dt * invmass * jX);
	LHS += identity;

	//Compute RHS  
	Eigen::VectorXd RHS(3*n); //right hand side hM^-1 * (fo + h * df/dX * v0) (from B&W) 		
	RHS = dt_ * invmass * (phi_zero + dt_ * jX * v_zero);

	time = omp_get_wtime() - time; 

	if(printThings) {
		cout << "setup LHS took : " << time << endl;
		
		cout << " LHS norm = " << LHS.norm() << " RHS norm = " << RHS.norm() << endl; 
	}

	time = omp_get_wtime();

	if(iterative) {	
		double time2 = omp_get_wtime(); 
		solver_it.analyzePattern(LHS);
		time2 = omp_get_wtime() - time2; 
		solver_it.factorize(LHS);
	}

	time = omp_get_wtime() - time; 
	if(printThings)	cout << "Total Factorize LHS took : " << time << endl; 


	time = omp_get_wtime(); 

	//SOLVE Step
	if(iterative) dv = solver_it.solveWithGuess(RHS,dv);

	time = omp_get_wtime() - time; 
	if(printThings)	cout << "Solve step took : " << time << endl; 

	int iters = 0;
	if(iterative) {
		iters = solver_it.iterations(); 
		if(printThings) cout << "solved in " << iters << " iterations and error = " << solver_it.error() << endl; 
	}

	//Update positions and velocities (new only) 
	Eigen::VectorXd dX(3*n);
	for(int i=0;i<n;i++) {
		Point* pi = points[i]; 
		Array dvi = Array(dv[3*i],dv[3*i+1],dv[3*i+2]);
		pi->vnew = pi->v + dvi; 		
		Array dXi = dt_ * pi->vnew;
		pi->Xnew = pi->X + dXi;
		pi->applyBC(dXi); 
		dX[3*i] = dXi[0]; 
		dX[3*i+1] = dXi[1]; 
		dX[3*i+2] = dXi[2]; 
	}
	
	updateInternalForces();
	Eigen::VectorXd phi_new(3*n);
	for(int i=0; i<n; i++) {
		Point* pi = points[i]; 
		Array forces = pi->phi + pi->fex; 	
		phi_new[3*i] = forces[0]; 
		phi_new[3*i+1] = forces[1]; 
		phi_new[3*i+2] = forces[2]; 
	}	
	
	//We use the taylor series of order zero to compute the error
	Eigen::VectorXd OdX_approx(3*n);
	OdX_approx = phi_new - (phi_zero + jX * dX); //this should be a O(dX) 
		
	double err0 = (OdX_approx).norm(); 

	for(int i=0; i<n; i++) {
		Point* pi = points[i]; 
		Array locerr(OdX_approx(3*i),OdX_approx(3*i+1),OdX_approx(3*i+2));
		pi->err = locerr.norm(); 
	}

	err = tol * err0;

	if(printThings) {
		cout << "dX norm = " << dX.norm() << endl;  
		cout << "phi_zero norm = " << phi_zero.norm() << endl; 
		cout << "computed error = " << err << endl; 
		cout << "err0 = " << err0 << endl; 
	}

	//Computation of dtnew
	double facold = facold_BW;
	double betadt = stab_factor; 
	double expo1 = 1./8. - betadt * 0.2;
	double fac11 = pow(err,expo1);	
	double fac1 = 0.333;
	double fac2 = 6.0;
	double facc1 = 1.0 / fac1;
	double facc2 = 1.0 / fac2;
	double safe = 0.9;
	double fac = fac11 / pow(facold,betadt);
	fac = std::max(facc2,std::min(facc1,fac/safe));
	double dtnew = dt_ / fac;

    if(use_variable_timestep) dt_ = dtnew; 
	
	if(!use_variable_timestep || ((err <= 1.))) {
		//Update positions and velocities
		for(int i=0;i<n;i++) {
			Point* pi = points[i]; 
	
			Array dummy(0.,0.,0.); 
			pi->applyBC(dummy); 	

			pi->v = pi->vnew; 
			pi->X = pi->Xnew;
		}
			
		//Update res and apply BC
		res = 0.;
		for(int i=0;i<n;i++) {
			Point* pi = points[i]; 	
			
			Array res_ = pi->fex + pi->phi; 
			pi->applyBC(res_);
			pi->res = res_.norm();
			res += pi->res;
		}		
	}
	else { //STEP IS REJECTED
		dt_ = 0.9 * dt_;
	//	tol = 1.1 * tol;
		bw_init = true;
		if(iters==max_iters) dt_ = 0.5 * dt_; 
		accepted = false;
		//Do NOT Update positions and velocities
		for(int i=0;i<n;i++) {
			Point* pi = points[i]; 

			pi->vnew = pi->v; 		 
			pi->Xnew = pi->X;	
		}
	}
	return accepted; 
}

void StructuralSolver::computeJacobian(spm& jX, spm& jV){
	if(delingette) {
		computeJacobianTRQS(jX,jV);
	}
	else if(argyris) {
		computeJacobianArgyris(jX,jV); 
	}
}

void StructuralSolver::computeJacobianTRQS(spm& jX,spm& jV) {
	int n = points.size(); 
	std::vector<T> jAl; 

	//Compute the position jacobian dF/dX
	for(size_t i=0; i<facets.size(); i++) {
		facets[i]->computeJacobianTRQS(jAl);
	}

	jX.setFromTriplets(jAl.begin(),jAl.end());	
	jX.makeCompressed();
	jAl.clear(); 	

	//Compute the velocity jacobian dF/dV
	if(beta_damp < 0.) {
        for(int i=0; i<n; i++) {
            Point* p = points[i];
                int id = p->index;	
            double coef = p->beta;
                if(!p->dirich){	
                jAl.push_back(T(3*id,3*id,coef));
                jAl.push_back(T(3*id+2,3*id+2,coef));
                if(!p->dirichY){
                    jAl.push_back(T(3*id+1,3*id+1,coef));
                }	
            }
        }
	}

	//visco elastic damping
	if(ksi < 0.) {
        for(size_t i=0; i<edges.size(); i++) {
            Edge* edge = edges[i]; 
            
            int id1 = edge->points.first->index; 
            int id2 = edge->points.second->index; 
            
            double coef = edge->damp_coef2; 
                
            if(!edge->points.first->dirich && !edge->points.second->dirich){	
                jAl.push_back(T(3*id1,3*id2,-coef));
                jAl.push_back(T(3*id1+2,3*id2+2,-coef));
        
                jAl.push_back(T(3*id2,3*id1,-coef));
                jAl.push_back(T(3*id2+2,3*id1+2,-coef));
                
                //diagonal terms
                jAl.push_back(T(3*id1,3*id1,coef));
                jAl.push_back(T(3*id1+2,3*id1+2,coef));
        
                jAl.push_back(T(3*id2,3*id2,coef));
                jAl.push_back(T(3*id2+2,3*id2+2,coef));
                
                if(!edge->points.first->dirichY && !edge->points.second->dirichY){
                    jAl.push_back(T(3*id1+1,3*id2+1,-coef));
                    jAl.push_back(T(3*id2+1,3*id1+1,-coef));
                    jAl.push_back(T(3*id1+1,3*id1+1,coef));
                    jAl.push_back(T(3*id2+1,3*id2+1,coef));
                }	
            }
        }
	}

	jV.setFromTriplets(jAl.begin(),jAl.end());	
	jV.makeCompressed();
}

void StructuralSolver::computeJacobianArgyris(spm& jX, spm& jV) {
	//-- Computes all the terms on the points (outside lines)
	preComputeJacobianArgyris(); 

	//-- Dispatch them into an Eigen SparseMatrix
	vector<T> jAl; 
	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i]; 
		int id1 = p->index; 

		//Diag terms. 
		Matrix3 diag = p->jacDiagTermArgy; 
		for(int k=0; k<3; k++) {
			for(int l=0; l<3; l++) {
				jAl.push_back(T(3*id1+k,3*id1+l,diag(k,l)));		
			}
		}

		//Non Diag terms. 
		for(size_t j=0; j< p->jacNonDiagTermArgy.size(); j++) {
			pair<Matrix3,int> elt = p->jacNonDiagTermArgy[j]; 
			Matrix3 ndiag = elt.first;
			int id2 = elt.second; 

			for(int k=0; k<3; k++) {
				for(int l=0; l<3; l++) {
					jAl.push_back(T(3*id1+k,3*id2+l,ndiag(k,l)));		
				}
			}
		}
	}

	jX.setFromTriplets(jAl.begin(),jAl.end());	
	jX.makeCompressed();

	jAl.clear(); 	

	//Compute the velocity jacobian dF/dV
	if(beta_damp < 0.) {
		for(size_t i=0; i<points.size(); i++) {
			Point* p = points[i];
			int id = p->index;	
			double coef = p->beta;
			if(!p->dirich) {	
				jAl.push_back(T(3*id,3*id,coef));
				jAl.push_back(T(3*id+2,3*id+2,coef));
				if(!p->dirichY){
					jAl.push_back(T(3*id+1,3*id+1,coef));
				}	
			}
		}
	}
	jV.setFromTriplets(jAl.begin(),jAl.end());	
	jV.makeCompressed();
}

void StructuralSolver::preComputeJacobianArgyris() {
    size_t i=0; 
    #pragma omp parallel private(i) 
    {	
        //Compute facet contributions
        #pragma omp for	
        for(i=0; i<facets.size(); i++) {
            facets[i]->computeJacobianArgyris(); 	
        }
    
        #pragma omp barrier 
    
        //Reduce diagonal terms on all points
        #pragma omp for	
        for(i=0; i<points.size(); i++) {
            points[i]->reduceJacobiansArgyris(); 	
        }
        
        #pragma omp barrier 

        //Reduce non diagonal terms on edges and
        //distribute them on the points. 
        #pragma omp for	
        for(i=0; i<edges.size(); i++) {
            edges[i]->reduceJacobiansArgyris(); 
        } 
    }  
}

void StructuralSolver::updateInternalForces() {
	if(delingette) {
		updateInternalForcesDelingette();
	}	
	else if(argyris) {
		updateInternalForcesArgyris();
	}
}

void StructuralSolver::updateInternalForcesDelingette() {
	size_t i = 0;
	#pragma omp parallel private(i)
	{
		#pragma omp for //Update edge lengths
		for(i=0; i<edges.size(); i++) {
			Edge* edge = edges[i]; 
			edge->c = 0.; //Reset shear coeff
			edge->computeLength(); 
		}
		Array zero = Array(0.,0.,0.);
		#pragma omp for 
		for(i=0; i<points.size(); i++) {
			points[i]->phi = zero; 
		}

		#pragma omp for  
		for(i=0; i<facets.size(); i++) {

			Facet* facet = facets[i]; //deformed triangle Tq	

			double dl1 = facet->edge1->getdl();
			double dl2 = facet->edge2->getdl();
			double dl3 = facet->edge3->getdl();
			
			//version with cvectors. Needs a reduction after this loop!
			facet->edge1->cvector[facet->edge1_cvector_index] = facet->c2*dl3 + facet->c3*dl2; 
			facet->edge2->cvector[facet->edge2_cvector_index] = facet->c3*dl1 + facet->c1*dl3;
			facet->edge3->cvector[facet->edge3_cvector_index] = facet->c1*dl2 + facet->c2*dl1;	
		}
	
		#pragma omp barrier

		#pragma omp for 
		for(i=0; i<edges.size(); i++) {
			Edge* edge = edges[i];
			
			//Reduction, only necessary if have the cvector version.
			//Comment this if we are in the lock version. 
			edge->reduceCvector(); 

			double invl = 1./edge->l; 
			
			//Simple relative velocity ksi damping term.
			Array damp_vec =  - damp_relax_param * edge->damp_coef2 * (edge->dV );//- dV.dot(normal) * normal);
			
			//Axial sigma damping term. 
			double damp = edge->damp_coef * edge->dV.dot(edge->dX) * invl;

			Array phi = (edge->k*edge->dl  + edge->c - damp) * invl * edge->dX + damp_vec; 

			//Version with phi vectors. 
			edge->points.first->phivector[edge->p1_phivector_index] = phi; 
			edge->points.second->phivector[edge->p2_phivector_index] = - phi;
		}

		#pragma omp barrier
		
		#pragma omp for 
		for(i=0; i<points.size(); i++) {
			Point* p = points[i]; 			
			//Reduction, only necessary in phiVector version. 	
			p->reducePhiVector(); 
			p->phi += damp_relax_param * p->beta * p->vnew;
		
		}
	
		#pragma omp barrier
	} //parallel section
}

void StructuralSolver::updateInternalForcesArgyris(){
	    size_t i = 0;

    #pragma omp parallel private(i)
    {
        //Update edge lengths
        #pragma omp for 	
        for(i=0;i<edges.size();i++){
            edges[i]->computeLength();
        }	
        #pragma omp for	
        for(i=0; i<facets.size(); i++) {
            Facet* facet = facets[i]; 
            facet->computeLength(); 
        }
    
        Array zero = Array(0.,0.,0.);
        #pragma omp for 
        for(i=0; i<points.size(); i++) {
            Point * p = points[i];
            p->phi = zero; 
        }
    
        /*---------------------Update forces with updateInternalForcesArgyris from facet class------------------*/
        #pragma omp barrier

        #pragma omp for 
        for(i=0; i<facets.size(); i++) {
            facets[i]->updateInternalForcesArgyris();
        }
        
        #pragma omp barrier

        /*---------------------Adding the damping-----------*/
        
        #pragma omp for
        for(i=0; i<points.size(); i++) {
            Point * p = points[i];
            p->reducePhiVector();
            
            p->phi += damp_relax_param * p->beta * p->vnew;

        }
        }//Jai fermé la section parallele ici pour voir 

        for(i=0; i<edges.size(); i++) {
            Edge* edge = edges[i];
            
            double invl = 1./edge->l; 	
            
            //Simple relative velocity ksi damping term. 
            Array damp_vec =  - damp_relax_param * edge->damp_coef2 * (edge->dV );//- dV.dot(normal) * normal);
            
            //	cout << "damp_vec " << damp_vec << endl;	
            //Axial sigma damping term. 
            double damp = edge->damp_coef * edge->dV.dot(edge->dX) * invl;
        
            //	cout << "damp : " << damp << endl;
            Array damp_to_add = ( - damp * invl * edge->dX+  damp_vec); 
            //	cout << "damp_to add / phi : " << damp_to_add << endl;
            //Version with phi vectors. 
            edge->points.first->phi += damp_to_add; 
            edge->points.second->phi += - damp_to_add;
    }
}

void StructuralSolver::computeMassesAndAreas() {
	
	//Set all masses and areas to zero
	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i]; 
		p->area = 0.;
		p->mass = 0.;	
	}
	
	double totalSurface = 0.0;
	double totalMass = 0.;
	for(size_t i=0; i<facets.size(); i++) {
		Point* p1 = facets[i]->points[0];
		Point* p2 = facets[i]->points[1];
		Point* p3 = facets[i]->points[2];
		double area = Utils::getArea(p1->X,p2->X,p3->X); 
		double mass = area * surfaceDensity; 
		totalSurface += area;	
		totalMass += mass; 
		area = 1.* area/3.; 
		mass = 1.* mass/3.;
		p1->area += area;
	    p1->mass += mass;
		p2->area += area;
	    p2->mass += mass;
		p3->area += area;
	    p3->mass += mass;	
	}
	
	int n = points.size(); 
	double avg_mass = 1. * totalSurface * surfaceDensity / n;

	vector<Array> null_mass_points; 

	for(size_t i=0; i<points.size(); i++) {
		Point* p = points[i];
		p->invmass = 1./p->mass;
		if(p->mass<1.e-14) {
			p->mass = avg_mass; 
			p->invmass = 1./avg_mass;
			null_mass_points.push_back(p->X);
		}
	}

	if(null_mass_points.size()>0) {	
		std::cout << "ERROR ! Found " << null_mass_points.size() << " points with null mass, printed for debug." << endl; 
		string debug_file = "null_mass_points.vtk"; 	
		Utils::printPoints(debug_file,null_mass_points); 
	}
	std::cout <<  " avg mass = " << avg_mass << std::endl; 
	cout << "Total surface = " << totalSurface << endl; 
	cout << "Total mass = " << totalMass << endl; 
}

void StructuralSolver::computeStresses() {
	max_strain = 0.;
	total_area = 0.;

	int index = 0;
	for(size_t i=0; i<facets.size(); i++) {
		Facet* facet = facets[i]; 
		
		Point* p1 = facet->points[0];
		Point* p2 = facet->points[1]; 
		Point* p3 = facet->points[2]; 
		
		double area = Utils::getArea(p1->X,p2->X,p3->X);	
		
		facets[index]->svm_max = 0;
	
		if (argyris){
			facets[i]->computeCartesianStressArgyris(); 
		}else{
			facets[i]->computeCartesianStress();
		}
		//Computing total area
		total_area += area;
	}
	for(size_t i=0; i<edges.size(); i++){
		double eps = edges[i]->getStrain();
		if (eps > max_strain){ 
			max_strain = eps;
		}
	}
}