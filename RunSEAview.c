#include <vector>
#include "plothelper.h"
#include "SEAview_lite_epem.h"
#include "SEAobject.h"

void RunSEAview(){

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
    TTreeFormula* f_true_opang = new TTreeFormula("TrueOpeningAngle","acos((mctruth_daughters_px[0]*mctruth_daughters_px[1]+mctruth_daughters_py[0]*mctruth_daughters_py[1]+mctruth_daughters_pz[0]*mctruth_daughters_pz[1])/(sqrt(mctruth_daughters_px[0]*mctruth_daughters_px[0]+mctruth_daughters_py[0]*mctruth_daughters_py[0]+mctruth_daughters_pz[0]*mctruth_daughters_pz[0])*sqrt(mctruth_daughters_px[1]*mctruth_daughters_px[1]+mctruth_daughters_py[1]*mctruth_daughters_py[1]+mctruth_daughters_pz[1]*mctruth_daughters_pz[1])))*180.0/3.14109",v);
    
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
    TH2D* h = new TH2D("True:Reco Opening Angle Response", "True:Reco Opening Angle Response",45,0,45,45,0,45);

    //some configuration bits
    std::vector<int> cols = {kBlue-6, kMagenta+1};
    bool do_plot_2d = true;
    double radius = 8.0;

    int cnt = 0; 
    //Loop over all entries
    for(size_t i=0; i< v->GetEntries(); i++){

        v->GetEntry(i);
        if(i%500==0)std::cout<<i<<"/"<<v->GetEntries()<<std::endl;

        //Right now focus on 2 object collections
        if( (num_reco_showers+num_reco_tracks>2)) continue;

        //This is necessary to get around a bug in root to do with vectors in TTrees. Call GetNdata() before evaluating TTreeFormula. 
        for(auto &f: forms) f->GetNdata();


        double easy = std::max(f_true_ep_E->EvalInstance(), f_true_em_E->EvalInstance())/(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance());
        if(easy>0.9) continue;

        double e_tot = (f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance());
        if(e_tot<0.1) continue;


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
            std::cout<<"On traco number "<<t<<" which has "<<reco_track_spacepoint_x->at(t).size()<<" "<<reco_track_spacepoint_y->at(t).size()<<" "<<reco_track_spacepoint_z->at(t).size()<<" spacepoints"<<std::endl;
            objs.emplace_back(1,cols[objs.size()], reco_track_hit_wire->at(t), reco_track_hit_plane->at(t), reco_track_hit_tick->at(t), reco_track_hit_energy->at(t), reco_track_spacepoint_x->at(t), reco_track_spacepoint_y->at(t), reco_track_spacepoint_z->at(t));
            objs.back().print();
        }


        double reco_ang = 0;

        //Some info to plot on the SEAviewer
        std::vector<std::string> tags ={"True #theta e^{+}e^{-}","E_max/E_total","E_total", "Num Trk: ","Num Shr: "};
        std::vector<std::string> vals = {  to_string_prec(true_opang,1),to_string_prec( std::max(f_true_ep_E->EvalInstance(), f_true_em_E->EvalInstance())/(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance())  ,2),to_string_prec(f_true_em_E->EvalInstance()+f_true_ep_E->EvalInstance(),3),std::to_string(num_reco_tracks),std::to_string(num_reco_showers)};


        //And do the cal
        reco_ang = SEAviewer(objs, reco_vertex_3D, reco_vertex_2D, true_vertex, std::to_string(subrun_number)+"_"+std::to_string(event_number), tags, vals, radius,do_plot_2d); 

        //Fill output Th2D
        h->Fill( reco_ang,true_opang);
        std::cout<<"How did we do? True OpAng : "<<true_opang<<" Reco OpAng : "<<reco_ang<<std::endl;

        cnt++;
        if(cnt>100)break;
   //     break;        
    }

    //Normalize the TH2D into a respsonse matrix.
    std::vector<double> norms;
    for(int i=0; i<=h->GetNbinsX(); i++){
        norms.push_back(0.0);
        for(int j=0; j<=h->GetNbinsY(); j++){
            norms.back()+= h->GetBinContent(i,j);
        }
    }
    for(int i=0; i<=h->GetNbinsX(); i++){
        for(int j=0; j<=h->GetNbinsY(); j++){
            h->SetBinContent(i,j, h->GetBinContent(i,j)/norms[i]);
        }
    
    }
    
    TCanvas *ch = new TCanvas();
    ch->cd();
    h->Draw("colz");
    h->GetYaxis()->SetTitle("Reco e^{+}e^{-} Opening Angle [Deg]");
    h->GetXaxis()->SetTitle("True e^{+}e^{-} Opening Angle [Deg]");
    ch->SaveAs("Response.pdf","pdf");
return;

}
