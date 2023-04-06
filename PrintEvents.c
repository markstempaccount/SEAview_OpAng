#include <string>
#include <iostream>
#include <stdio.h>
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"

int PrintEvents(){
	std::string work_dir = "/uboone/data/users/ltong/";
	std::string data_dir = "PointClouds";
	gSystem -> cd(work_dir.c_str());
	gSystem -> mkdir(data_dir.c_str());
	gSystem -> cd(data_dir.c_str());
	TFile *fWC = new TFile("/pnfs/uboone/persistent/users/markross/IsoTropic_epem_samples/reprocessed_wc_pandora_epem_file_v2.root", "READ");
	TTree *v = (TTree*) fWC -> Get("newtree;1");
	int npan_sp, nwc_sp;
	v->SetBranchAddress("npan_sp", &npan_sp);
	v->SetBranchAddress("nwc_sp", &nwc_sp);
	int max_sp = 5000;
	double pan_sp_x[max_sp];
	double pan_sp_y[max_sp];
	double pan_sp_z[max_sp];
	double wc_sp_x[max_sp];
	double wc_sp_y[max_sp];
	double wc_sp_z[max_sp];
	v->SetBranchAddress("pan_sp_x",&pan_sp_x);
	v->SetBranchAddress("pan_sp_y",&pan_sp_y);
	v->SetBranchAddress("pan_sp_z",&pan_sp_z);
	v->SetBranchAddress("wc_sp_x",&wc_sp_x);
	v->SetBranchAddress("wc_sp_y",&wc_sp_y);
	v->SetBranchAddress("wc_sp_z",&wc_sp_z);
	Long64_t run_number=0;
        Long64_t subrun_number=0;
       	Long64_t event_number=0;
        int in_wirecell = 1;
        v->SetBranchAddress("run",&run_number);
        v->SetBranchAddress("subrun",&subrun_number);
        v->SetBranchAddress("event",&event_number);
	double true_opang;
        v -> SetBranchAddress("true_opening_angle", &true_opang);
	double reco_pan_vertex_x, reco_pan_vertex_y, reco_pan_vertex_z, reco_wc_vertex_x, reco_wc_vertex_y, reco_wc_vertex_z;
        v -> SetBranchAddress("reco_pan_vertex_x", &reco_pan_vertex_x);
        v -> SetBranchAddress("reco_pan_vertex_y", &reco_pan_vertex_y);
        v -> SetBranchAddress("reco_pan_vertex_z", &reco_pan_vertex_z);
	v -> SetBranchAddress("reco_wc_vertex_x", &reco_wc_vertex_x);
        v -> SetBranchAddress("reco_wc_vertex_y", &reco_wc_vertex_y);
        v -> SetBranchAddress("reco_wc_vertex_z", &reco_wc_vertex_z);
	double true_vertex_x, true_vertex_y, true_vertex_z;
        v -> SetBranchAddress("true_wc_vertex_x", &true_vertex_x);
        v -> SetBranchAddress("true_wc_vertex_y", &true_vertex_y);
        v -> SetBranchAddress("true_wc_vertex_z", &true_vertex_z);
	double true_ep_E, true_em_E;
        v -> SetBranchAddress("true_positron_energy", &true_ep_E);
        v -> SetBranchAddress("true_electron_energy", &true_em_E);
	for(int i=0; i< v->GetEntries(); i++){
		if(i%100 == 0)
			std::cout << "Begin Processing " << i << "th Event\n";
		v->GetEntry(i);
		if(npan_sp == 0 && nwc_sp == 0){
			std::cout << "Skipped Run " << run_number << "Subrun " << subrun_number << "Event " << event_number << "due to no spacepoints.\n";
			continue;
		}
		std::string filename = "PointCloud_" + std::to_string(run_number) +  "_" + std::to_string(subrun_number) + "_" + std::to_string(event_number) + ".dat";
		FILE *file;
		file = fopen(filename.c_str(), "w");
		fprintf(file, "Run: %lld\n", run_number);
		fprintf(file, "Subrun: %lld\n", subrun_number);
		fprintf(file, "Event: %lld\n", event_number);
		fprintf(file, "True Electron Energy (MeV): %lf\n", true_em_E);
		fprintf(file, "True Positron Energy (MeV): %lf\n", true_ep_E);
		fprintf(file, "True Opening Angle (degrees): %lf\n", true_opang);
		fprintf(file, "True Vertex (cm): %lf %lf %lf\n", true_vertex_x, true_vertex_y, true_vertex_z);
		fprintf(file, "Reconstructed Wirecell Vertex (cm): %lf %lf %lf\n", reco_wc_vertex_x, reco_wc_vertex_y, reco_wc_vertex_z);
		fprintf(file, "Reconstructed Pandora Vertex (cm): %lf %lf %lf\n\n", reco_pan_vertex_x, reco_pan_vertex_y, reco_pan_vertex_z);
		fprintf(file, "Number of Wirecell Spacepoints: %d\n", nwc_sp);
		fprintf(file, "Wirecell Spacepoints\n");
		for(int j = 0; j < nwc_sp; j++)
			fprintf(file, "%lf %lf %lf\n", wc_sp_x[j], wc_sp_y[j], wc_sp_z[j]);
		fprintf(file, "\n");
		fprintf(file, "Number of Pandora Spacepoints: %d\n", npan_sp);
		fprintf(file, "Pandora Spacepoints\n");
		for(int j = 0; j < npan_sp; j++)
			fprintf(file, "%lf %lf %lf\n", pan_sp_x[j], pan_sp_y[j], pan_sp_z[j]);
		fclose(file);
	}
	return 0;
}
