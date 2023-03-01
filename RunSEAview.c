//C++ includes
#include <vector>
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <map>

//Our local includes
#include "plothelper.h"
#include "SEAview_lite_epem.h"
#include "SEAobject.h"
#include "SetRealAspectRatio2.h"

//Other includes
#include "CLI11.hpp"

//Root includes
#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <Fit/Fitter.h>
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
#include "TRandom3.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"


int main(int argc, char **argv){

    CLI::App app{"Run SEAviewer for e+e- opening angle calculations"}; 

    // Define options
    double radius = 1.0; bool candles = true; bool subplots = false; bool graphEVD_SEAview = false; bool graphResponse = true; bool normalizeResponse = false; double dist_2_true_max = 1.0; double min_rtang_diff = 0; int maxcnt = 1e6; double easymax = 0.9; double e_totmin = 100; bool iterateRadius = true; double radiusInterval = 0.5; double maxradius = 20;

    //doubles
    app.add_option("-r,--radius", radius, "radius for search");
    app.add_option("-d,--dist", dist_2_true_max, "distance from true to reco");
    app.add_option("-m,--mindiff", min_rtang_diff, "Min reco-trie ang");
    app.add_option("-n,--maxcnt", maxcnt, "max count of TTree entries");
    app.add_option("-e,--easymax", easymax, "Energy Asymetry max");
    app.add_option("-t,--emin", e_totmin, "Total energy min");
    app.add_option("-i,--interval", radiusInterval,"radius to iterative over");
    app.add_option("-x,--maxradius", maxradius,"Max radius");

    //bools
    app.add_option("-c,--candles", candles, "Plot candles?");
    app.add_option("-s,--subplots",subplots, "Plot subplots?");
    app.add_option("-v,--evd",graphEVD_SEAview,"Plot EVDs?");
    app.add_option("-p,--response",graphResponse,"Plot EVDs?");
    app.add_option("-z,--normalize",normalizeResponse,"Normalize Response?");
    app.add_option("-a,--iterate",iterateRadius,"Iterate Radius?");

    CLI11_PARSE(app, argc, argv);

    printf("Argument input for this run -- radius : %f , dist : %f , mindiff : %f, maxcnt : %i, easymax : %f, emin : %f, interval : %f, maxrad : %f \n",radius,dist_2_true_max,min_rtang_diff,maxcnt,easymax,e_totmin,radiusInterval,maxradius);
    printf("Argument bools for this run -- candles : %i , subplots : %i , evd : %i, response : %i, normalize : %i, iterate: %i \n",candles,subplots,graphEVD_SEAview,graphResponse,normalizeResponse,iterateRadius);


    //If not iterating radius, set maxradius to current radius so main loop only runs once.
    if(!iterateRadius)
        maxradius = radius;
    if(radiusInterval <= 0){
        std::cout << "radiusInterval cannot be negative or 0!" << std::endl;
        return 0;
    }

    std::string sp_reco = "wc"; //wc or pan
    std::string vtx_reco = "pan";

    //Some vectors to store the 2D recob::Hits. The Wire number, the plane (0,1 or 2), the energy and the peak time tick
    std::vector<std::vector<int>> *reco_shower_hit_wire = new std::vector<std::vector<int>>(2, std::vector<int>(2, 1));
    std::vector<std::vector<int>> *reco_shower_hit_plane = new std::vector<std::vector<int>>(2, std::vector<int>(2, 1));
    std::vector<std::vector<double>> *reco_shower_hit_energy = new std::vector<std::vector<double>>(2, std::vector<double>(2, 1));
    std::vector<std::vector<double>> *reco_shower_hit_tick = new std::vector<std::vector<double>>(2, std::vector<double>(2, 1));
    std::vector<std::vector<int>> *reco_track_hit_wire = new std::vector<std::vector<int>>(2, std::vector<int>(2, 1));
    std::vector<std::vector<int>> *reco_track_hit_plane = new std::vector<std::vector<int>>(2, std::vector<int>(2, 1));
    std::vector<std::vector<double>> *reco_track_hit_energy =new std::vector<std::vector<double>>(2, std::vector<double>(2, 1));
    std::vector<std::vector<double>> *reco_track_hit_tick = new std::vector<std::vector<double>>(2, std::vector<double>(2, 1));

    //h6 are histograms of error vs. radius. Only used if iterating through radii but must be declared regardless to have sufficient scope.
    TH2D *h6absolute, *h6percent;
    if(iterateRadius){
        h6absolute = new TH2D("h6absolute", "Radius:Reco Opening Angle Error", (int) ((maxradius - radius)/radiusInterval) + 1, radius - radiusInterval/2, maxradius + radiusInterval/2, 225, -45, 180);
        h6percent = new TH2D("h6percent", "Radius:Reco Opening Angle Error", (int) ((maxradius - radius)/radiusInterval) + 1, radius - radiusInterval/2, maxradius + radiusInterval/2, 100, -100, 200);
    }

    while(radius <= maxradius){

        //base directory the input files are in
        std::string base_dir = "/pnfs/uboone/persistent/users/markross/Jan2022_gLEE_files/";
        //std::string base_dir = "/home/mark/work/uBooNE/EplusEmin_Retreat2022_Master/files/";

        //Parent directory of data director(y/ies)
        std::string work_dir;
        if(sp_reco == "wc")
            work_dir = "/uboone/app/users/ltong/eplus_eminus_studies_2023/SEAview_OpAng/wirecell";
        else if(sp_reco == "pan")
            work_dir = "/uboone/app/users/ltong/eplus_eminus_studies_2023/SEAview_OpAng";

        //Data directory. Keeps output from different radii separate, which is useful if iterating. Unfortunately, requires regenerating dictionaries for each radius.
        std::string response_dir = (std::to_string(radius));			

        //Remove trailing 0s from directory name
        response_dir.erase ( response_dir.find_last_not_of('0') + 1, std::string::npos );
        response_dir.erase ( response_dir.find_last_not_of('.') + 1, std::string::npos );
        response_dir = response_dir + "cm";

        //Make directory for current radius and cd to it.
        gSystem -> cd(work_dir.c_str());
        gSystem -> mkdir(response_dir.c_str());
        gSystem -> cd(response_dir.c_str());

        //Grab the TTrees associated with the gLEE Ntuples (automatically friends the necessary bits internally). Function in plothelper.h simple and quick. 
        TFile *fWC = new TFile("/pnfs/uboone/persistent/users/markross/IsoTropic_epem_samples/reprocessed_wc_pandora_epem_file_v2.root", "READ");
        TTree *v = (TTree*) fWC -> Get("newtree;1");

        //Create some dictionaryies so root can read vectors of vectors. Not strictly needed but safer. 
        //		gInterpreter->GenerateDictionary("std::vector<std::vector<int> >", "vector");
        //		gInterpreter->GenerateDictionary("std::vector<std::vector<double> >", "vector");

        //Don't bother showing things to screen, just pring PDF
        gROOT->SetBatch(kTRUE);

        //Some vectors to store spacepoints. Spacepoints are the 3D reconstructed objects corresponding to the "matched" 2D recob::Hits on each plane
        /*std::vector<std::vector<double>> *wc_sp_x =0;
          std::vector<std::vector<double>> *wc_sp_y =0;
          std::vector<std::vector<double>> *wc_sp_z =0;
          v->SetBranchAddress("wc_sp_x",&wc_sp_x);
          v->SetBranchAddress("wc_sp_y",&wc_sp_y);
          v->SetBranchAddress("wc_sp_z",&wc_sp_z);
          std::vector<std::vector<double>> *wc_sp_x =0;
          std::vector<std::vector<double>> *wc_sp_y =0;
          std::vector<std::vector<double>> *wc_sp_z =0;
          v->SetBranchAddress("wc_sp_x",&wc_sp_x);
          v->SetBranchAddress("wc_sp_y",&wc_sp_y);
          v->SetBranchAddress("wc_sp_z",&wc_sp_z);*/

        int num_sp;
        v->SetBranchAddress(("n"+sp_reco+"_sp").c_str(), &num_sp);
        int max_sp = 5000;
        double sp_x[max_sp];
        double sp_y[max_sp];
        double sp_z[max_sp];
        v->SetBranchAddress((sp_reco+"_sp_x").c_str(),&sp_x);
        v->SetBranchAddress((sp_reco+"_sp_y").c_str(),&sp_y);
        v->SetBranchAddress((sp_reco+"_sp_z").c_str(),&sp_z);


        //Some Other info for the event. How many tracks and showers pandora neutrino slice reconstructed, as well as Run:Subrun:Event number for ID
        int num_reco_showers = 2;
        int num_reco_tracks = 0;
        Long64_t run_number=0;
        Long64_t subrun_number=0;
        Long64_t event_number=0;
        int in_wirecell = 1;
        v->SetBranchAddress("run",&run_number);
        v->SetBranchAddress("subrun",&subrun_number);
        v->SetBranchAddress("event",&event_number);

        //If you want to make complex formula from many TTree branches, using SetBranchAddress becomes tedious. 
        //So we can also access more complex formula directly via TTreeFormula and not individual branches (often betteir)
        //E.g the true opening angle between the e+e- pair.
        //For these isotropic files, the mctruth_daughters_XX is a vector of length 2 with [0] and [1] corresponding to the e+e- pair. 
        double true_opang;
        v -> SetBranchAddress("true_opening_angle", &true_opang);

        //Vertex point of vertex in 3D
        double reco_vertex_x, reco_vertex_y, reco_vertex_z;
        v -> SetBranchAddress(("reco_"+vtx_reco+"_vertex_x").c_str(), &reco_vertex_x);
        v -> SetBranchAddress(("reco_"+vtx_reco+"_vertex_y").c_str(), &reco_vertex_y);
        v -> SetBranchAddress(("reco_"+vtx_reco+"_vertex_z").c_str(), &reco_vertex_z);

        //Start point of vertex in each 2D wire plane (wire,tick)
        std::vector<double> reco_vertex_2D = {1,1,1,1};

        /*Direction of shower in 3D //Not currenrly used really
          TTreeFormula* f_reco_shower_direction_x = new TTreeFormula("reco_shower_direction_x","reco_shower_dirx[0]",v);
          TTreeFormula* f_reco_shower_direction_y = new TTreeFormula("reco_shower_direction_y","reco_shower_diry[0]",v);
          TTreeFormula* f_reco_shower_direction_z = new TTreeFormula("reco_shower_direction_z","reco_shower_dirz[0]",v);*/

        //Some truth level information. Little trickier for TextGen as we dont save a "nu vertex" to space charge correct, so we rely on the sim info. 
        double true_vertex_x, true_vertex_y, true_vertex_z;
        v -> SetBranchAddress("true_wc_vertex_x", &true_vertex_x);
        v -> SetBranchAddress("true_wc_vertex_y", &true_vertex_y);
        v -> SetBranchAddress("true_wc_vertex_z", &true_vertex_z);
        //TTreeFormula* f_true_vertex_x = new TTreeFormula("true_vertex_x","mctruth_daughters_startx[0]",v);
        //TTreeFormula* f_true_vertex_y = new TTreeFormula("true_vertex_y","mctruth_daughters_starty[0]",v);
        //TTreeFormula* f_true_vertex_z = new TTreeFormula("true_vertex_z","mctruth_daughters_startz[0]",v);


        /*True Direction of epluz and eminus (unit vectors)
          TTreeFormula* f_true_ep_dir_x = new TTreeFormula("true_eplus_dir_x","mctruth_daughters_px[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
          TTreeFormula* f_true_ep_dir_y = new TTreeFormula("true_eplus_dir_y","mctruth_daughters_py[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
          TTreeFormula* f_true_ep_dir_z = new TTreeFormula("true_eplus_dir_z","mctruth_daughters_pz[0]/sqrt(pow(mctruth_daughters_px[0],2)+pow(mctruth_daughters_py[0],2)+pow(mctruth_daughters_pz[0],2))",v);
          TTreeFormula* f_true_em_dir_x = new TTreeFormula("true_eminus_dir_x","mctruth_daughters_px[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);
          TTreeFormula* f_true_em_dir_y = new TTreeFormula("true_eminus_dir_y","mctruth_daughters_py[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);
          TTreeFormula* f_true_em_dir_z = new TTreeFormula("true_eminus_dir_z","mctruth_daughters_pz[1]/sqrt(pow(mctruth_daughters_px[1],2)+pow(mctruth_daughters_py[1],2)+pow(mctruth_daughters_pz[1],2))",v);*/

        //True eplus and einus energies
        double true_ep_E, true_em_E;
        v -> SetBranchAddress("true_positron_energy", &true_ep_E);
        v -> SetBranchAddress("true_electron_energy", &true_em_E);

        //This is just a safety issue to make a necessary bug fix. See a few lines down
        //std::vector<TTreeFormula*> forms = {f_true_opang,f_reco_vertex_x,f_reco_vertex_y,f_reco_vertex_z,f_reco_shower_direction_x,f_reco_shower_direction_y,f_reco_shower_direction_z,f_true_vertex_x,f_true_vertex_y,f_true_vertex_z, f_true_ep_dir_x,    f_true_ep_dir_y,    f_true_ep_dir_z,    f_true_em_dir_x,    f_true_em_dir_y,    f_true_em_dir_z,f_reco_vertex_wire_plane0,f_reco_vertex_wire_plane1,f_reco_vertex_wire_plane2,f_reco_vertex_tick, f_true_ep_E, f_true_em_E};


        //Make Histograms to show relation between reco opening angle error and various reco and true parameters.
        TH2D *h1 = new TH2D("h1", "True:Reco Opening Angle Response",45,0,45,90,0,90);
        TH2D *h2absolute = new TH2D("h2absolute", "Num Tracks + Showers:Reco Opening Angle Error", 2, 0.5, 2.5, 135, -45, 90);
        TH2D *h3absolute = new TH2D("h3absolute", "E_max/E_total:Reco Opening Angle Error", (int) std::round((easymax - 0.50)/0.01), 0.50, easymax, 135, -45, 90);
        TH2D *h4absolute = new TH2D("h4absolute", "True:Reco Opening Angle Error", 45, 0, 45, 135, -45, 90);
        TH2D *h5absolute = new TH2D("h5absolute", "Reco Vertex Error:Reco Opening Angle Error", (int) std::round(dist_2_true_max/0.025), 0, dist_2_true_max, 135, -45, 90);
        TH2D *h2percent = new TH2D("h2percent", "Num Tracks + Showers:Reco Opening Angle Error", 2, 0.5, 2.5, 100, -100, 200);
        TH2D *h3percent = new TH2D("h3percent", "E_max/E_total:Reco Opening Angle Error", (int) std::round((easymax - 0.50)/0.01), 0.50, easymax, 100, -100, 200);
        TH2D *h4percent = new TH2D("h4percent", "True:Reco Opening Angle Error", 45, 0, 45, 100, -100, 200);
        TH2D *h5percent = new TH2D("h5percent", "Reco Vertex Error:Reco Opening Angle Error", (int) std::round(dist_2_true_max/0.025), 0, dist_2_true_max, 100, -100, 200);

        //some configuration bits
        std::vector<int> cols = {kBlue-6, kMagenta+1};
        bool do_plot_2d = false;

        int cnt = 0; 
        //Loop over all entries
        std::cout<<"We have :"<<v->GetEntries()<<" in the file."<<std::endl;
        for(size_t i=0; i< v->GetEntries(); i++){

            std::cout<<"***************Start**************** "<<i<<std::endl;
            std::cout<<i<<"/"<<v->GetEntries()<<std::endl;

            v->GetEntry(i);

            std::vector<double> *reco_spacepoint_x = new std::vector<double>(sp_x, sp_x + num_sp);
            std::vector<double> *reco_spacepoint_y = new std::vector<double>(sp_y, sp_y + num_sp);
            std::vector<double> *reco_spacepoint_z = new std::vector<double>(sp_z, sp_z + num_sp);

            bool shall_we = true;

            if(!in_wirecell || reco_spacepoint_x->size()==0 ) shall_we = false;
            //Right now focus on 2 object collections
            if((num_reco_showers == 0) || (num_reco_showers+num_reco_tracks>2)) shall_we = false;

            //This is necessary to get around a bug in root to do with vectors in TTrees. Call GetNdata() before evaluating TTreeFormula. 
            //for(auto &f: forms) f->GetNdata();

            //e_assymetry shouldn't be too large, or program can't distinguish that there's 2 showers.
            double easy = std::max(true_ep_E, true_em_E)/(true_em_E+true_ep_E);
            if(easy>easymax) shall_we = false;

            //E_total shouldn't be too small, or there won't be proper trails.
            double e_tot = (true_em_E+true_ep_E);
            if(e_tot<e_totmin) shall_we = false;

            //some variables for easy use, can always access the formula directily too (x,y,z)
            std::vector<double> reco_vertex_3D = {reco_vertex_x,reco_vertex_y,reco_vertex_z};
            std::vector<double> true_vertex =  {true_vertex_x,true_vertex_y,true_vertex_z};

            //Distance between Pandora/Wirecell reco vertex and true vertex. If too large, can't properly reconstruct tracks and showers.
            double dist_2_true = sqrt( pow(reco_vertex_3D[0]-true_vertex[0],2)+  pow(reco_vertex_3D[1]-true_vertex[1],2) + pow(reco_vertex_3D[2]-true_vertex[2],2) );
            std::cout<<"Dist: "<<dist_2_true<<", Reco: "<<reco_vertex_3D[0]<<" "<<reco_vertex_3D[1]<<" "<<reco_vertex_3D[2]<<" True: "<<true_vertex[0]<<" "<<true_vertex[1]<<" "<<true_vertex[2]<<std::endl;
            if(dist_2_true > dist_2_true_max) shall_we = false;


            if(shall_we){
                //Build up a vector of "Objects", tracks showers and unassociated hits
                std::vector<SEAobject> objs;

                /*if(false){
                //std::cout<<"Starting to set up "<<num_reco_showers<<" ( "<<wc_sp_z->size()<<" ) showers and "<<num_reco_tracks<<" ( "<<wc_sp_z->size()<<" ) tracks . in event "<<i<<std::endl;
                //showers
                for(int s=0; s< num_reco_showers; s++){
                //std::cout<<"On Shower number "<<s<<" which has "<<wc_sp_y->at(s).size()<<" spacepoints"<<std::endl;
                objs.emplace_back(0, cols[objs.size()], reco_shower_hit_wire->at(s), reco_shower_hit_plane->at(s), reco_shower_hit_tick->at(s), reco_shower_hit_energy->at(s), wc_sp_x->at(s), wc_sp_y->at(s), wc_sp_z->at(s));
                objs.back().print();
                std::cout<<"Hmm"<<std::endl;
                }

                //tracks
                for(int t=0; t< num_reco_tracks; t++){
                std::cout<<"On track number "<<t<<" which has "<<wc_sp_x->at(t).size()<<" "<<wc_sp_y->at(t).size()<<" "<<wc_sp_z->at(t).size()<<" spacepoints"<<std::endl;
                objs.emplace_back(1,cols[objs.size()], reco_track_hit_wire->at(t), reco_track_hit_plane->at(t), reco_track_hit_tick->at(t), reco_track_hit_energy->at(t), wc_sp_x->at(t), wc_sp_y->at(t), wc_sp_z->at(t));
                objs.back().print();
                }
                }else{
                */
                if(reco_track_hit_wire->size()!=0){
                    objs.emplace_back(1, cols[0], reco_track_hit_wire->at(0), reco_track_hit_plane->at(0), reco_track_hit_tick->at(0), reco_track_hit_energy->at(0), *reco_spacepoint_x, *reco_spacepoint_y, *reco_spacepoint_z);
                }else if(reco_track_hit_wire->size()!=0){
                    objs.emplace_back(1, cols[0], reco_shower_hit_wire->at(0), reco_shower_hit_plane->at(0), reco_shower_hit_tick->at(0), reco_shower_hit_energy->at(0), *reco_spacepoint_x, *reco_spacepoint_y, *reco_spacepoint_z);
                }else {
                    std::vector<int> tmpi = {0};
                    std::vector<double> tmpd ={0};
                    objs.emplace_back(1, cols[0], tmpi,tmpi,tmpd,tmpd, *reco_spacepoint_x, *reco_spacepoint_y, *reco_spacepoint_z);

                }


                //}


                //double reco_ang = 0;

                //Some info to plot on the SEAviewer
                std::vector<std::string> tags ={"True #theta e^{+}e^{-}","E_max/E_total","E_total", "Num Trk: ","Num Shr: "};
                std::vector<std::string> vals = {  to_string_prec(true_opang,1),to_string_prec(easy,2),to_string_prec(e_tot,3),std::to_string(num_reco_tracks),std::to_string(num_reco_showers)};

                //And do the cal does calculation and plots
                //SEAviewer reco_obj (objs, reco_vertex_3D, reco_vertex_2D, true_vertex, sp_reco + "_sp_" + vtx_reco + "_vert_"+std::to_string(subrun_number)+"_"+std::to_string(event_number), tags, vals, radius); 
                SEAviewer reco_obj (objs, true_vertex, reco_vertex_2D, true_vertex, "fixed_plane_fit_"+std::to_string(subrun_number)+"_"+std::to_string(event_number), tags, vals, radius); 
                reco_obj.reco_ang_calc();

                //Whether or not to graph individual events.
                if(graphEVD_SEAview)
                    reco_obj.plotter();


                //Want to look at failure modes; set min_rtang_diff to 0 (default) if wishing to view all events.
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
                    if(iterateRadius){
                        //std::cout<<1.1<<std::endl;
                        h6absolute -> Fill(radius, absoluteerror);
                        //std::cout<<1.2<<std::endl;
                        h6percent -> Fill(radius, percenterror);
                        //std::cout<<1.3<<std::endl;
                    }
                }
                std::cout<<"How did we do? True OpAng : "<<true_opang<<" Reco OpAng : "<<reco_obj.reco_ang<<std::endl;
            }

            delete reco_spacepoint_x;
            delete reco_spacepoint_y;
            delete reco_spacepoint_z;

            cnt++;
            if(cnt>maxcnt) break; //Number of pdfs to save
        }

        //Normalize (or not) the TH2D into a respsonse matrix.
        if(graphResponse){
            std::string fResponsename = sp_reco + "Response" + "trueVertex";
            TH2D* h[9] = {h1, h2absolute, h2percent, h3absolute, h3percent, h4absolute, h4percent, h5absolute, h5percent};
            std::string xlabel[5] = {"True e^{+}e^{-} Opening Angle [Deg]", "Number of Tracks + Showers", "E_max/E_total", "True e^{+}e^{-} Opening Angle [Deg]", "Vertex Error [cm]"}; //Is Vertex error actually in cm?
            TFile *fResponse = new TFile((fResponsename + ".root").c_str(), "RECREATE");
            fResponse -> cd();
            gDirectory -> mkdir(response_dir.c_str());
            fResponse -> cd(response_dir.c_str());
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
                //h[k] -> GetYaxis() -> SetRangeUser(firstbinylowedge, lastbinyhighedge);
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
                    SetRealAspectRatio2(ch1, 1);
                    SetRealAspectRatio2(ch1, 2);
                    ch1 ->Write();
                    ch1 ->SaveAs((fResponsename + ".pdf(").c_str(), "pdf");
                    ch1 ->SaveAs("ch1.C");
                }
                else if(k == 8){
                    ch3 ->Write();
                    ch3 -> SaveAs((fResponsename + ".pdf)").c_str(), "pdf");
                    ch3 ->SaveAs("ch3.C");
                }
                else if(!subplots){
                    ch3 ->Write();
                    ch3 -> SaveAs((fResponsename + ".pdf").c_str(), "pdf");
                    ch3 ->SaveAs("ch3.C");
                }
                else if(k == 7){
                    ch2 -> Write();
                    ch3 -> SaveAs((fResponsename + ".pdf").c_str(), "pdf");
                    ch2 -> SaveAs("ch2.C");
                }
                h[k] ->SaveAs((fResponsename + std::to_string(k) + ".C").c_str());
                fResponse -> Close();
                //delete ch1;
                //delete ch2;
                //delete ch3;
            }}
        gSystem -> cd(work_dir.c_str());
        radius += radiusInterval;
        fWC->Close();
        //delete h1;
        //delete h2absolute;
        //delete h2percent;
        //delete h3absolute;
        //delete h3percent;
        //delete h4absolute;
        //delete h4percent;
        //delete h5absolute;
        //delete h5percent;
    }
    if(iterateRadius){
        TH2D* h[2] = {h6absolute, h6percent};	
        TCanvas *ch = new TCanvas("ch2", "ch2", 3000, 2400);
        for(int k = 0; k < 2; k++){
            std::vector<double> norms;
            for(int i=0; i<=h[k] ->GetNbinsX(); i++){
                norms.push_back(0.0);
                for(int j=0; j<=h[k] ->GetNbinsY(); j++){
                    norms.back()+= h[k] ->GetBinContent(i,j);
                }}
            for(int i=0; i<=h[k] ->GetNbinsX(); i++){
                for(int j=0; j<=h[k] ->GetNbinsY(); j++){
                    h[k] ->SetBinContent(i,j, h[k] ->GetBinContent(i,j)/norms[i]);
                }}
            ch -> cd();
            int firstbiny = h[k] -> FindFirstBinAbove(0, 2);
            double firstbinylowedge = h[k] -> GetYaxis() -> GetBinLowEdge(firstbiny);
            int lastbiny = h[k] -> FindLastBinAbove(0, 2);
            double lastbinyhighedge = h[k] -> GetYaxis() -> GetBinLowEdge(lastbiny + 1);
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
            if(k==0)
                h[k] -> GetYaxis() -> SetTitle("Reco e^{+}e^{-} Opening Angle Error [Deg]");
            else
                h[k] -> GetYaxis() -> SetTitle("Reco e^{+}e^{-} Opening Angle Error [%]");
            h[k] -> GetXaxis() ->SetTitle("Radius of Reco Circle (cm)");
            if(k == 0)
                ch -> SaveAs((sp_reco + "Radius" + "trueVertex.pdf(").c_str(), "pdf");
            else
                ch -> SaveAs((sp_reco + "Radius" + "trueVertex.pdf)").c_str(), "pdf");
        }}

    delete reco_shower_hit_wire;
    delete reco_shower_hit_plane;
    delete reco_shower_hit_energy ;
    delete reco_shower_hit_tick ;
    delete reco_track_hit_wire ;
    delete reco_track_hit_plane;
    delete reco_track_hit_energy;
    delete reco_track_hit_tick ;



    return 0;

}
