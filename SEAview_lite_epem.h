#ifndef SEAVIEWLITE_H
#define SEAVIEWLITE_H


#include "SEAobject.h"
#include "RecoOpAng1.h"
#include "TH1D.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"
#include "TLegend.h"
#include "TPolyLine.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <map>
#include "TRandom3.h"
#include "TLine.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <algorithm>

#include <string>


//Create a SEAviewer class, to show 2D and 3D hits and spacepoints visually
class SEAviewer
{

	public:
	//Variables passed to object
	std::vector<SEAobject> &objs; 
	std::vector<double> &reco_vertex_3D; 
	std::vector<double> &reco_vertex_2D; 
	std::vector<double> &true_vertex; 
	std::string tag; 
	std::vector<std::string> tags; 
	std::vector<std::string> vals;
	double radius;

	//Variables related to graphing. Not ottally sure if I need to initalize these here or not	
	std::vector<double> all_fit_points_x;
	std::vector<double> all_fit_points_y;
	std::vector<double> all_fit_points_z;
	std::vector<double> all_fit_weights;
	std::vector<TGraph> out_graphs2D; //to save output left/right pts
	std::vector<double> left_fit;//the parameters of the fitted lines, for visualization only
	std::vector<double> right_fit;

	//Output
	double reco_ang = -9999;


	//Constructor to initalize variables
	SEAviewer(  std::vector<SEAobject> &objs1, 
                std::vector<double> &reco_vertex_3D1, 
                std::vector<double> &reco_vertex_2D1,   
                std::vector<double> &true_vertex1, 
                std::string tag1, std::vector<std::string> tags1,
                std::vector<std::string>vals1, 
                double radius1) :objs(objs1), reco_vertex_3D(reco_vertex_3D1), reco_vertex_2D(reco_vertex_2D1), true_vertex(true_vertex1), tag(tag1), tags(tags1), vals(vals1), radius(radius1) {}


	int  reco_ang_calc();

	void plotter();


};
#endif
