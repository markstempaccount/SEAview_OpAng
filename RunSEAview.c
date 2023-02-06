#include <vector>
#include "plothelper.h"
#include "SEAview_lite_epem.h"
#include "SEAobject.h"
#include "SetRealAspectRatio.h"

void RunSEAview(double radius = 8.0, bool candles = true, bool subplots = false, bool graphEVD_SEAview = false, bool graphResponse = true, bool normalizeResponse = false, double dist_2_true_max = 1.0, double min_rtang_diff = 0, int maxcnt = 1e6, double easymax = 0.9, double e_totmin = 0.1, bool iterateRadius = false, int loop = 1){

	//base directory the input files are in
	std::string base_dir = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/";
	//std::string base_dir = "/home/mark/work/uBooNE/EplusEmin_Retreat2022_Master/files/";

	//Grab the TTrees associated with the gLEE Ntuples (automatically friends the necessary bits internally). Function in plothelper.h simple and quick. 
	TTree *v = (TTree*)loadgLEE(base_dir+"vertex_Isotropic_EpEm_Batch1_v50.5_SP.root", "singlephotonana");

	//Create some dictionaryies so root can read vectors of vectors. Not strickly needed but safer. 
	gInterpreter->GenerateDictionary("std::vector<std::vector<int> >", "vector");
	gInterpreter->GenerateDictionary("std::vector<std::vector<double> >", "vector");

	//Don't bother showing things to screen, just pring PDF
	gROOT->SetBatch(kTRUE);

	//Some vectors to store spacepoints. Spacepoints are the 3D reconstricted objects corresponding to the "matched" 2D recob::Hits on each plane
	std::vector<std::vector<double>> *reco_shower_spacepoint_x =0;
	std::vector<std::vector<double>> *reco_shower_spacepoint_y =0;
	std::vector<std::vector<double>> *reco_shower_spacepoint_z =0;
	v->SetBranchAddress("reco_shower_spacepoint_x",&reco_shower_spacepoint_x);
	v->SetBranchAddress("reco_shower_spacepoint_y",&reco_shower_spacepoint_y);
	v->SetBranchAddress("reco_shower_spacepoint_z",&reco_shower_spacepoint_z);
	std::vector<std::vector<double>> *reco_track_spacepoint_x =0;
	std::vector<std::vector<double>> *reco_track_spacepoint_y =0;
	std::vector<std::vector<double>> *reco_track_spacepoint_z =0;
	v->SetBranchAddress("reco_track_spacepoint_x",&reco_track_spacepoint_x);
	v->SetBranchAddress("reco_track_spacepoint_y",&reco_track_spacepoint_y);
	v->SetBranchAddress("reco_track_spacepoint_z",&reco_track_spacepoint_z);

	//Some vectors to store the 2D recob::Hits. The Wire number, the plane (0,1 or 2), the energy and the peak time tick
	std::vector<std::vector<int>> *reco_shower_hit_wire = 0;
	std::vector<std::vector<int>> *reco_shower_hit_plane = 0;
	std::vector<std::vector<double>> *reco_shower_hit_energy = 0;
	std::vector<std::vector<double>> *reco_shower_hit_tick = 0;
	v->SetBranchAddress("reco_shower_hit_wire",&reco_shower_hit_wire);
	v->SetBranchAddress("reco_shower_hit_tick",&reco_shower_hit_tick);
	v->SetBranchAddress("reco_shower_hit_plane",&reco_shower_hit_plane);
	v->SetBranchAddress("reco_shower_hit_energy",&reco_shower_hit_energy);
	std::vector<std::vector<int>> *reco_track_hit_wire = 0;
	std::vector<std::vector<int>> *reco_track_hit_plane = 0;
	std::vector<std::vector<double>> *reco_track_hit_energy = 0;
	std::vector<std::vector<double>> *reco_track_hit_tick = 0;
	v->SetBranchAddress("reco_track_hit_wire",&reco_track_hit_wire);
	v->SetBranchAddress("reco_track_hit_tick",&reco_track_hit_tick);
	v->SetBranchAddress("reco_track_hit_plane",&reco_track_hit_plane);
	v->SetBranchAddress("reco_track_hit_energy",&reco_track_hit_energy);


	//Some Other info for the event. How many tracks and showers pandora neutrino slice reconstructed, as well as Run:Subrun:Event number for ID
	int num_reco_showers = 0;
	int num_reco_tracks = 0;
	int run_number = 0;
	int subrun_number = 0;
	int event_number = 0;
	v->SetBranchAddress("reco_asso_showers",&num_reco_showers);
	v->SetBranchAddress("reco_asso_tracks",&num_reco_tracks);
	v->SetBranchAddress("run_number",&run_number);
	v->SetBranchAddress("subrun_number",&subrun_number);
	v->SetBranchAddress("event_number",&event_number);

	//If you want to make complex formula from many TTree branches, using SetBranchAddress becomes tedious. 
	//So we can also access more complex formula directly via TTreeFormula and not individual branches (often betteir)
	//E.g the true opening angle between the e+e- pair.
	//For these isotropic files, the mctruth_daughters_XX is a vector of length 2 with [0] and [1] corresponding to the e+e- pair. 
	TTreeFormula* f_true_opang = new TTreeFormula("TrueOpeningAngle","acos((mctruth_daughters_px[0]*mctruth_daughters_px[1]+mctruth_daughters_py[0]*mctruth_daughters_py[1]+mctruth_daughters_pz[0]*mctruth_daughters_pz[1])/(sqrt(mctruth_daughters_px[0]*mctruth_daughters_px[0]+mctruth_daughters_py[0]*mctruth_daughters_py[0]+mctruth_daughters_pz[0]*mctruth_daughters_pz[0])*sqrt(mctruth_daughters_px[1]*mctruth_daughters_px[1]+mctruth_daughters_py[1]*mctruth_daughters_py[1]+mctruth_daughters_pz[1]*mctruth_daughters_pz[1])))*180.0/TMath::Pi()",v);

	//Vertex point of vertex in 3D
	TTreeFormula* f_reco_vertex_x = new TTreeFormula("reco_vertex_x","reco_vertex_x",v);
	TTreeFormula* f_reco_vertex_y = new TTreeFormula("reco_vertex_y","reco_vertex_y",v);
	TTreeFormula* f_reco_vertex_z = new TTreeFormula("reco_vertex_z","reco_vertex_z",v);

	//Start point of vertex in each 2D wire plane (wire,tick)
	TTreeFormula* f_reco_vertex_wire_plane0 = new TTreeFormula("reco_vertex_wire_plane0","reco_vertex_wire_p0",v);
	TTreeFormula* f_reco_vertex_wire_plane1 = new TTreeFormula("reco_vertex_wire_plane1","reco_vertex_wire_p1",v);
	TTreeFormula* f_reco_vertex_wire_plane2 = new TTreeFormula("reco_vertex_wire_plane2","reco_vertex_wire_p2",v);
	TTreeFormula* f_reco_vertex_tick = new TTreeFormula("reco_vertex_tick","reco_vertex_tick",v);

	//Direction of shower in 3D //Not currenrly used really
	TTreeFormula* f_reco_shower_direction_x = new TTreeFormula("reco_shower_direction_x","reco_shower_dirx[0]",v);
	TTreeFormula* f_reco_shower_direction_y = new TTreeFormula("reco_shower_direction_y","reco_shower_diry[0]",v);
	TTreeFormula* f_reco_shower_direction_z = new TTreeFormula("reco_shower_direction_z","reco_shower_dirz[0]",v);

	//Some truth level information. Little trickier for TextGen as we dont save a "nu vertex" to space charge correct, so we rely on the sim info. 
	TTreeFormula* f_true_vertex_x = new TTreeFormula("true_vertex_x","Alt$(sim_shower_start_x[0],Min$(sim_track_startx))",v);
	TTreeFormula* f_true_vertex_y = new TTreeFormula("true_vertex_y","Alt$(sim_shower_start_y[0],Min$(sim_track_starty))",v);
	TTreeFormula* f_true_vertex_z = new TTreeFormula("true_vertex_z","Alt$(sim_shower_start_z[0],Min$(sim_track_startz))",v);

	//True Direction of eplus and eminus (unit vectors)
	TTreeFormula* f_true_ep_dir_x = new TTreeFormula("true_eplus_dir_x","mctruth_daughters_px[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
	TTreeFormula* f_true_ep_dir_y = new TTreeFormula("true_eplus_dir_y","mctruth_daughters_py[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
	TTreeFormula* f_true_ep_dir_z = new TTreeFormula("true_eplus_dir_z","mctruth_daughters_pz[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
	TTreeFormula* f_true_em_dir_x = new TTreeFormula("true_eminus_dir_x","mctruth_daughters_px[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);
	TTreeFormula* f_true_em_dir_y = new TTreeFormula("true_eminus_dir_y","mctruth_daughters_py[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);
	TTreeFormula* f_true_em_dir_z = new TTreeFormula("true_eminus_dir_z","mctruth_daughters_pz[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);

	//True eplus and einus energies
	TTreeFormula* f_true_ep_E = new TTreeFormula("true_eplus_E","mctruth_daughters_E[0]",v);
	TTreeFormula* f_true_em_E = new TTreeFormula("true_eminus_E","mctruth_daughters_E[1]",v);

	//This is just a safety issue to make a necessary bug fix. See a few lines down
	std::vector<TTreeFormula*> forms = {f_true_opang,f_reco_vertex_x,f_reco_vertex_y,f_reco_vertex_z,f_reco_shower_direction_x,f_reco_shower_direction_y,f_reco_shower_direction_z,f_true_vertex_x,f_true_vertex_y,f_true_vertex_z, f_true_ep_dir_x,    f_true_ep_dir_y,    f_true_ep_dir_z,    f_true_em_dir_x,    f_true_em_dir_y,    f_true_em_dir_z,f_reco_vertex_wire_plane0,f_reco_vertex_wire_plane1,f_reco_vertex_wire_plane2,f_reco_vertex_tick, f_true_ep_E, f_true_em_E};


	//Make a Histogram to save the output 2D response
	TH2D *h1 = new TH2D("h1", "True:Reco Opening Angle Response",45,0,45,180,0,180);
	TH2D *h2absolute = new TH2D("h2absolute", "Num Tracks + Showers:Reco Opening Angle Error", 2, 1, 3, 225, -45, 180);
	TH2D *h3absolute = new TH2D("h3absolute", "E_max/E_total:Reco Opening Angle Error", (int) std::round((easymax - 0.50)/0.01), 0.50, easymax, 225, -45, 180);
	TH2D *h4absolute = new TH2D("h4absolute", "True:Reco Opening Angle Error", 45, 0, 45, 225, -45, 180);
	TH2D *h5absolute = new TH2D("h5absolute", "Reco Vertex Error:Reco Opening Angle Error", 40, 0, dist_2_true_max, 225, -45, 180);
	TH2D *h2percent = new TH2D("h2percent", "Num Tracks + Showers:Reco Opening Angle Error", 2, 1, 3, 100, -100, 200);
	TH2D *h3percent = new TH2D("h3percent", "E_max/E_total:Reco Opening Angle Error", (int) std::round((easymax - 0.50)/0.01), 0.50, easymax, 100, -100, 200);
	TH2D *h4percent = new TH2D("h4percent", "True:Reco Opening Angle Error", 45, 0, 45, 100, -100, 200);
	TH2D *h5percent = new TH2D("h5percent", "Reco Vertex Error:Reco Opening Angle Error", 40, 0, dist_2_true_max, 100, -100, 200);

//	TH2D *h6 = new TH2D("Radius:Reco Opening Angle Error", 

	//some configuration bits
	std::vector<int> cols = {kBlue-6, kMagenta+1};
	bool do_plot_2d = true;

	int cnt = 0; 
	//Loop over all entries
	std::cout<<"We have :"<<v->GetEntries()<<" in the file."<<std::endl;
	for(size_t i=0; i< v->GetEntries(); i++){

		std::cout<<"***************Start****************"<<std::endl;

		v->GetEntry(i);
		if(i%500==0)std::cout<<i<<"/"<<v->GetEntries()<<std::endl;

		//Right now focus on 2 object collections
		if((num_reco_showers+num_reco_tracks==0) || (num_reco_showers+num_reco_tracks>2)) continue;

		//This is necessary to get around a bug in root to do with vectors in TTrees. Call GetNdata() before evaluating TTreeFormula. 
		for(auto &f: forms) f->GetNdata();


		double easy = std::max(f_true_ep_E->EvalInstance(), f_true_em_E->EvalInstance())/(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance());
		if(easy>easymax) continue;

		double e_tot = (f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance());
		if(e_tot<e_totmin) continue;

		//Grab the instances for this event for some information
		double true_opang = f_true_opang->EvalInstance();

		//some variables for easy use, can always access the formula directily too (x,y,z)
		std::vector<double> reco_vertex_3D = {f_reco_vertex_x->EvalInstance(),f_reco_vertex_y->EvalInstance(),f_reco_vertex_z->EvalInstance() };
		std::vector<double> reco_vertex_2D = {f_reco_vertex_wire_plane0->EvalInstance(),f_reco_vertex_wire_plane1->EvalInstance(),f_reco_vertex_wire_plane2->EvalInstance() ,f_reco_vertex_tick->EvalInstance()};
		std::vector<double> true_vertex =  {f_true_vertex_x->EvalInstance(),f_true_vertex_y->EvalInstance(),f_true_vertex_z->EvalInstance() };


		double dist_2_true = sqrt( pow(reco_vertex_3D[0]-true_vertex[0],2)+  pow(reco_vertex_3D[1]-true_vertex[1],2) + pow(reco_vertex_3D[2]-true_vertex[2],2) );
		std::cout<<"Dist: "<<dist_2_true<<", Reco: "<<reco_vertex_3D[0]<<" "<<reco_vertex_3D[1]<<" "<<reco_vertex_3D[2]<<" True: "<<true_vertex[0]<<" "<<true_vertex[1]<<" "<<true_vertex[2]<<std::endl;
		if(dist_2_true > 1) continue;

		//Build up a vector of "Objects", tracks showers and unassociated hits
		std::vector<SEAobject> objs;

		std::cout<<"Starting to set up "<<num_reco_showers<<" ( "<<reco_shower_spacepoint_z->size()<<" ) showers and "<<num_reco_tracks<<" ( "<<reco_track_spacepoint_z->size()<<" ) tracks . in event "<<i<<std::endl;
		//showers
		for(int s=0; s< num_reco_showers; s++){
			std::cout<<"On Shower number "<<s<<" which has "<<reco_shower_spacepoint_y->at(s).size()<<" spacepoints"<<std::endl;
			objs.emplace_back(0, cols[objs.size()], reco_shower_hit_wire->at(s), reco_shower_hit_plane->at(s), reco_shower_hit_tick->at(s), reco_shower_hit_energy->at(s), reco_shower_spacepoint_x->at(s), reco_shower_spacepoint_y->at(s), reco_shower_spacepoint_z->at(s));
			objs.back().print();
			std::cout<<"Hmm"<<std::endl;
		}

		//tracks
		for(int t=0; t< num_reco_tracks; t++){
			std::cout<<"On track number "<<t<<" which has "<<reco_track_spacepoint_x->at(t).size()<<" "<<reco_track_spacepoint_y->at(t).size()<<" "<<reco_track_spacepoint_z->at(t).size()<<" spacepoints"<<std::endl;
			objs.emplace_back(1,cols[objs.size()], reco_track_hit_wire->at(t), reco_track_hit_plane->at(t), reco_track_hit_tick->at(t), reco_track_hit_energy->at(t), reco_track_spacepoint_x->at(t), reco_track_spacepoint_y->at(t), reco_track_spacepoint_z->at(t));
			objs.back().print();
		}


		//double reco_ang = 0;

		//Some info to plot on the SEAviewer
		std::vector<std::string> tags ={"True #theta e^{+}e^{-}","E_max/E_total","E_total", "Num Trk: ","Num Shr: "};
		std::vector<std::string> vals = {  to_string_prec(true_opang,1),to_string_prec( std::max(f_true_ep_E->EvalInstance(), f_true_em_E->EvalInstance())/(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance())  ,2),to_string_prec(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance(),3),std::to_string(num_reco_tracks),std::to_string(num_reco_showers)};

		//And do the cal does calculation and plots
		SEAviewer reco_obj (objs, reco_vertex_3D, reco_vertex_2D, true_vertex, "bad_vertex_"+std::to_string(subrun_number)+"_"+std::to_string(event_number), tags, vals, radius); 

		reco_obj.reco_ang_calc();
		if(graphEVD_SEAview)
			reco_obj.plotter();


		//Want to look at failure modes
		float rtang_diff = fabs(reco_obj.reco_ang - true_opang);
		if(rtang_diff < min_rtang_diff) continue;
		std::cout<<"*******************************************Difference between True and Reco Ang "<<rtang_diff<<std::endl;	

		//Fill output Th2D
		if(graphResponse){
			h1 -> Fill(true_opang, reco_obj.reco_ang);
			double absoluteerror = reco_obj.reco_ang - true_opang;
			double percenterror = absoluteerror/true_opang*100;
			h2absolute -> Fill(num_reco_showers+num_reco_tracks, absoluteerror);
			h3absolute -> Fill(easy, absoluteerror);
			h4absolute -> Fill(true_opang, absoluteerror);
			h5absolute -> Fill(dist_2_true, absoluteerror);
			h2percent -> Fill(num_reco_showers+num_reco_tracks, percenterror);
			h3percent -> Fill(easy, percenterror);
			h4percent -> Fill(true_opang, percenterror);
			h5percent -> Fill(dist_2_true, percenterror);
		}
		std::cout<<"How did we do? True OpAng : "<<true_opang<<" Reco OpAng : "<<reco_obj.reco_ang<<std::endl;

		cnt++;
		if(cnt>maxcnt) break; //Number of pdfs to save
	}

	//Normalize the TH2D into a respsonse matrix.
	if(graphResponse){
		TH2D* h[9] = {h1, h2absolute, h2percent, h3absolute, h3percent, h4absolute, h4percent, h5absolute, h5percent};
		std::string xlabel[5] = {"True e^{+}e^{-} Opening Angle [Deg]", "Number of Tracks + Showers", "E_max/E_total", "True e^{+}e^{-} Opening Angle [Deg]", "Vertex Error [cm]"}; //Is Vertex error actually in cm?
		TCanvas *ch1 = new TCanvas();
		TCanvas *ch2 = new TCanvas("ch2", "ch2", 3000, 2400);
		TCanvas *ch3 = new TCanvas("ch3", "ch3", 3000, 2400);
		if(subplots){
			ch2 -> Divide(2, 2);
			ch3 -> Divide(2, 2);
		}
		for(int k = 0; k < 9; k++){
			if(normalizeResponse){
				std::vector<double> norms;
				for(int i=0; i<=h[k] ->GetNbinsX(); i++){
					norms.push_back(0.0);
					for(int j=0; j<=h[k] ->GetNbinsY(); j++){
						norms.back()+= h[k] ->GetBinContent(i,j);
					}
				}
				for(int i=0; i<=h[k] ->GetNbinsX(); i++){
					for(int j=0; j<=h[k] ->GetNbinsY(); j++){
						h[k] ->SetBinContent(i,j, h[k] ->GetBinContent(i,j)/norms[i]);
					}

				}}

			ch1 -> cd();
			if(k != 0){
				if(subplots){
					if(k%2)
						ch2 -> cd((int) k/2.0 + 0.5);
					else
						ch3 -> cd((int) k/2.0);
				}
				else
					ch3 -> cd();
			}
			double firstbinxlowedge = h[k] -> GetXaxis() -> GetBinLowEdge(1);
			int lastbinx = h[k] -> FindLastBinAbove(0, 1);
			double lastbinxhighedge = h[k] -> GetXaxis() -> GetBinLowEdge(lastbinx + 1);
			int firstbiny = h[k] -> FindFirstBinAbove(0, 2);
			double firstbinylowedge = h[k] -> GetYaxis() -> GetBinLowEdge(firstbiny);
			int lastbiny = h[k] -> FindLastBinAbove(0, 2);
			double lastbinyhighedge = h[k] -> GetYaxis() -> GetBinLowEdge(lastbiny + 1);
			h[k] -> GetXaxis() -> SetRangeUser(firstbinxlowedge, lastbinxhighedge);
			h[k] -> GetYaxis() -> SetRangeUser(firstbinylowedge, lastbinyhighedge);
			h[k] -> SetFillStyle(0);
			h[k] -> SetLineColor(2);
			h[k] -> SetLineWidth(1);
			if(candles){
				h[k] ->Draw("colz candlex(00000311)");
				//h[k] ->Draw("colz");
				//TProfile *p = h[k] -> ProfileX("p", 1, -1, "s");
				//p -> Draw("same");
			}
			else
				h[k] -> Draw("colz");
			if(k == 0)
				h[k] -> GetYaxis()->SetTitle("Reco e^{+}e^{-} Opening Angle [Deg]");
			else if(k%2)
				h[k] -> GetYaxis() -> SetTitle("Reco e^{+}e^{-} Opening Angle Error [Deg]");
			else
				h[k] -> GetYaxis() -> SetTitle("Reco e^{+}e^{-} Opening Angle Error [%]");
			h[k] ->GetXaxis()->SetTitle(xlabel[(int) (k/2.0 + 0.5)].c_str());
			if(k == 0){
				SetRealAspectRatio(ch1, 1);
				SetRealAspectRatio(ch1, 2);
				ch1 ->SaveAs("Response.pdf(","pdf");
				ch1 ->SaveAs("ch1.C");
			}
			else if(k == 8){
				ch3 ->SaveAs("Response.pdf)","pdf");
				ch3 ->SaveAs("ch3.C");
			}
			else if(!subplots){
				ch3 ->SaveAs("Response.pdf", "pdf");
				ch3 ->SaveAs("ch3.C");
			}
			else if(k == 7){
				ch2 -> SaveAs("Response.pdf", "pdf");
				ch2 -> SaveAs("ch2.C");
			}
			//ch -> SaveAs(("Response" + std::to_string(k) + ".pdf").c_str(), "pdf");
			h[k] ->SaveAs(("Response" + std::to_string(k) + ".C").c_str());
	}}
	return;

}
