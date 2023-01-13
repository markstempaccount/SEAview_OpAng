#include "SEAobject.h"
#include "RecoOpAng1.h"

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


	//*********Reco calc function here**********
	void  reco_ang_calc(){ 
	//first get things set up to calculate reco_ang
	//Will fill these with the fittable points.
   
    //Loop over all objects and only select those within Xcm of reco vertex
    std::cout<<"SEAviewer::Begininig to calculate "<<radius<<" cm hits from all objs"<<std::endl;
    for(auto &obj: objs){
        for(int j=0; j<obj.f_num_sp; j++){
            double dist_2_vert = sqrt( pow(reco_vertex_3D[0]-obj.f_sp_x[j],2)+  pow(reco_vertex_3D[1]-obj.f_sp_y[j],2) + pow(reco_vertex_3D[2]-obj.f_sp_z[j],2) );
            if(dist_2_vert < radius){
                all_fit_points_x.push_back(obj.f_sp_x[j]);
                all_fit_points_y.push_back(obj.f_sp_y[j]);
                all_fit_points_z.push_back(obj.f_sp_z[j]);
            }
        }
    }

    //Our outputs will be the angle calculated
        //We need minimum 2 points for this algorithm!
    if(all_fit_points_x.size()>1){
        reco_ang =   recoOpAng1(reco_vertex_3D,   all_fit_points_x,    all_fit_points_y,    all_fit_points_z,out_graphs2D,left_fit,right_fit);
    }else{
        std::cout<<"WARNING: only "<< all_fit_points_x.size()<<" points withint "<<radius<<" cm so cant calc angle, defaulting to -9999"<<std::endl;
    }


    //Get things set up to print to PDF if wanted, otherwise return the single number in degrees
    return;
	}//End of reco calc




    //********************************** everything below here is plotting only *********************************//

	void plotter(){
    std::string print_name = "EVD_SEAview_"+tag;
    TCanvas * can=new TCanvas(print_name.c_str(),print_name.c_str(),3000,2400);
    can->Divide(4,3,0,0.1);

    double plot_point_size=1.0;        
    double buffer = 0.05;

    int n_objs = objs.size();
    double real_tick_min = 99999;
    double real_tick_max = -99999;
    for(auto &o: objs){
        real_tick_min = std::min(real_tick_min, o.f_min_tick);
        real_tick_max = std::max(real_tick_max, o.f_max_tick);
    }

    std::cout<<"Max and Min Ticks "<<real_tick_max<<" "<<real_tick_min<<std::endl;


    //******************************* First plot "Vertex" ***************************************
    std::vector<TGraph> vertex_graphs;
    std::vector<TGraph> vec_all_graphs;
    for(int i=0; i< 3; i++){
        vertex_graphs.emplace_back(1,&(reco_vertex_2D)[i], &(reco_vertex_2D)[3]);
    }

    std::vector<double> real_wire_min(3,99999); //real x axis edges for 3 planes
    std::vector<double> real_wire_max(3,-99999);

    for(int i=0; i<3; i++){
        TPad * pader = (TPad*)can->cd(i+1);
        if(i==0) pader->SetLeftMargin(0.1);

        //only show area surrounding the vertex up to std::min(plot_distance, distance_bw_vertex_channel_min/max)
        for(auto &o: objs){
            real_wire_min[i] = std::min(real_wire_min[i], o.f_min_wires[i]); 
            real_wire_max[i] = std::max(real_wire_max[i], o.f_max_wires[i]); 
        }

        /////////////////////////////////////////////////////////////////

        //        real_energy_min[i] =(double)*min_element(all_energies[i].begin(), all_energies[i].end()); ; 
        //        real_energy_max[i] =(double)*max_element(all_energies[i].begin(), all_energies[i].end()); ; 

        std::cout<<"SEAview:: "<<i<<" "<<real_wire_min[i]<<" - "<<real_wire_max[i]<<std::endl;
        //    std::cout<<"SEAview:: "<<i<<" EN: "<<real_energy_min[i]<<" - "<<real_energy_max[i]<<std::endl;
        //    std::cout<<"SEAview:: "<<vertex_graphs[i].GetN()<<" "<<vec_all_graphs[i].GetN()<<std::endl;

        vertex_graphs[i].SetMarkerStyle(29);
        vertex_graphs[i].SetMarkerSize(3);
        vertex_graphs[i].SetMarkerColor(kMagenta-3);
        vertex_graphs[i].GetYaxis()->SetRangeUser((1.0-buffer)*real_tick_min,(1.0+buffer)*real_tick_max);
        vertex_graphs[i].GetXaxis()->SetLimits((1.0-buffer)*real_wire_min[i], (1.0+buffer)*real_wire_max[i]);
        vertex_graphs[i].SetTitle(("Plane " +std::to_string(i)).c_str());
        vertex_graphs[i].GetYaxis()->SetTitle("Peak Hit Time Tick");
        vertex_graphs[i].GetXaxis()->SetTitle( ("Wire Number Plane " +std::to_string(i)).c_str());
        vertex_graphs[i].Draw("ap");

        if(i>0){
            vertex_graphs[i].GetYaxis()->SetLabelOffset(999);
            vertex_graphs[i].GetYaxis()->SetLabelSize(0);
        }
    }

    /********************************* All Hits in all OBJ's****************************/


    for(int i=0; i<3; i++){
        can->cd(i+1);
        int o = 0;
        for(auto &obj: objs){

            if(obj.f_graph_2D[i].GetN()>0){//need a check in case this track has no hits on this plane.

                obj.f_graph_2D[i].Draw("p same"); 
                obj.f_graph_2D[i].SetMarkerColor(obj.f_col);
                //obj.f_graph_2D[i].SetFillColor(kWhite);
                obj.f_graph_2D[i].SetMarkerStyle(20);
                obj.f_graph_2D[i].SetMarkerSize(plot_point_size*0.9);
            }
            o++;
        }
    }

    //******************************** DeadWireRegions********************************************
    //*** For now ignore, might copy in later ****
    //******************************** DeadWireRegions********************************************
    /* for(size_t i=0; i< m_bad_channel_list.size(); i++){
       int badchan = m_bad_channel_list[i].first;                                       
       int ok = m_bad_channel_list[i].second;       

       if(ok>1)continue;
       auto hs = geom->ChannelToWire(badchan); //type of hs: vector containing the ID of all the connected wires

       int thisp = (int)hs[0].Plane;
       double bc = hs[0].Wire;

       if(real_wire_min[thisp] < bc && bc < real_wire_max[thisp] ){
    //if(chan_min[thisp]-chan_shift < bc && bc < chan_max[thisp]+chan_shift ){
    can->cd(thisp+1);
    TLine *l = new TLine(bc, real_tick_min, bc, real_tick_max);
    //TLine *l = new TLine(bc,tick_min-tick_shift,bc,tick_max+tick_shift);
    l->SetLineColor(kGray+1);
    l->Draw("same");
    //can->cd(thisp+5);// Guanqun: how many values can plane ID take?
    //l->Draw("same");
    //can->cd(thisp+9);
    //l->Draw("same");
    }
    }
    */


    // The 3D stuff
    //**************************** INFO ***************************/

    double X_max = -9999; 
    double Y_max = -9999;
    double Z_max = -9999;

    double X_min = 9999; 
    double Y_min = 9999;
    double Z_min = 9999;

    for(auto &obj:objs){
        X_max = std::max(X_max, obj.f_max_sp[0]);
        Y_max = std::max(Y_max, obj.f_max_sp[1]);
        Z_max = std::max(Z_max, obj.f_max_sp[2]);

        X_min = std::min(X_min, obj.f_min_sp[0]);
        Y_min = std::min(Y_min, obj.f_min_sp[1]);
        Z_min = std::min(Z_min, obj.f_min_sp[2]);
    }

    double X_width = fabs(X_max-X_min);
    double Y_width = fabs(Y_max-Y_min);
    double Z_width = fabs(Z_max-Z_min);

    double X_mid = X_min+X_width/2.0;
    double Y_mid = Y_min+Y_width/2.0;
    double Z_mid = Z_min+Z_width/2.0;


    double ymodl = -1;
    double ymodu = 1;
    if(Y_min <0 ) ymodl = 1;
    if(Y_max <0)  ymodu = -1;

    TPad * paderXY = (TPad*)can->cd(5);
    paderXY->SetLeftMargin(0.1);
    TGraph  evdXY =   objs[0].f_graph_3D[0]; 
    TGraph * evdRecoStartPointXY = new TGraph(1,&(reco_vertex_3D)[0], &(reco_vertex_3D)[1]);
    TGraph * evdTrueStartPointXY = new TGraph(1,&(true_vertex)[0], &(true_vertex)[1]);
    evdXY.SetMarkerColor(objs[0].f_col);
    evdXY.SetMarkerStyle(20);
    evdXY.SetMarkerSize(1);
    evdXY.SetTitle("");
    evdXY.GetYaxis()->SetTitle(" Y [cm]");
    evdXY.GetXaxis()->SetTitle(" X [cm]");
    evdXY.GetYaxis()->SetRangeUser((1.0+ymodl*buffer)*Y_min,(1.0+ymodu*buffer)*Y_max);
    evdXY.GetXaxis()->SetLimits((1.0-buffer)*X_min, (1.0+buffer)*X_max);
    evdXY.Draw("ap");
    for(int i=1; i<objs.size(); i++){
        objs[i].f_graph_3D[0].SetMarkerStyle(20);
        objs[i].f_graph_3D[0].SetMarkerSize(1);
        objs[i].f_graph_3D[0].SetMarkerColor(objs[i].f_col);
        objs[i].f_graph_3D[0].Draw("p same");
    }

    evdRecoStartPointXY->Draw("same p");
    evdRecoStartPointXY->SetMarkerStyle(29);
    evdRecoStartPointXY->SetMarkerSize(3);
    evdRecoStartPointXY->SetMarkerColor(kRed);
    evdTrueStartPointXY->Draw("same p");
    evdTrueStartPointXY->SetMarkerStyle(29);
    evdTrueStartPointXY->SetMarkerSize(2);
    evdTrueStartPointXY->SetMarkerColor(kCyan);

    TEllipse *elXY3 = new TEllipse( reco_vertex_3D[0],reco_vertex_3D[1],3,3);
    elXY3->SetFillStyle(0);
    elXY3->SetLineColor(kRed-6);
    elXY3->SetLineWidth(1);
    elXY3->Draw();
    TEllipse *elXY5 = new TEllipse( reco_vertex_3D[0],reco_vertex_3D[1],radius,radius);
    elXY5->SetFillStyle(0);
    elXY5->SetLineColor(kMagenta);
    elXY5->SetLineWidth(1);
    elXY5->Draw();

    TPad * paderZY = (TPad*)can->cd(6);
    TGraph  evdZY =  objs[0].f_graph_3D[1];//new TGraph(reco_shower_sp_z.size(), &(reco_shower_sp_z)[0], &(reco_shower_sp_y)[0]);
    TGraph * evdRecoStartPointZY = new TGraph(1,&(reco_vertex_3D)[2], &(reco_vertex_3D)[1]);
    TGraph * evdTrueStartPointZY = new TGraph(1,&(true_vertex)[2], &(true_vertex)[1]);
    evdZY.SetMarkerStyle(20);
    evdZY.SetMarkerSize(1);
    evdZY.SetMarkerColor(objs[0].f_col);
    evdZY.SetTitle("");
    evdZY.GetYaxis()->SetTitle(" Y [cm]");
    evdZY.GetXaxis()->SetTitle(" Z [cm]");
    evdZY.GetYaxis()->SetRangeUser((1.0+ymodl*buffer)*Y_min,(1.0+ymodu*buffer)*Y_max);
    evdZY.GetXaxis()->SetLimits((1.0-buffer)*Z_min, (1.0+buffer)*Z_max);
    evdZY.Draw("ap");

    for(int i=1; i<objs.size(); i++){
        objs[i].f_graph_3D[1].SetMarkerStyle(20);
        objs[i].f_graph_3D[1].SetMarkerSize(1);
        objs[i].f_graph_3D[1].SetMarkerColor(objs[i].f_col);
        objs[i].f_graph_3D[1].Draw("p same");
    }


    evdRecoStartPointZY->Draw("same p");
    evdRecoStartPointZY->SetMarkerStyle(29);
    evdRecoStartPointZY->SetMarkerSize(3);
    evdRecoStartPointZY->SetMarkerColor(kRed);
    evdTrueStartPointZY->Draw("same p");
    evdTrueStartPointZY->SetMarkerStyle(29);
    evdTrueStartPointZY->SetMarkerSize(2);
    evdTrueStartPointZY->SetMarkerColor(kCyan);

    TEllipse *elZY3 = new TEllipse( reco_vertex_3D[2],reco_vertex_3D[1],3,3);
    elZY3->SetFillStyle(0);
    elZY3->SetLineColor(kRed-6);
    elZY3->SetLineWidth(1);
    elZY3->Draw();
    TEllipse *elZY5 = new TEllipse( reco_vertex_3D[2],reco_vertex_3D[1],radius,radius);
    elZY5->SetFillStyle(0);
    elZY5->SetLineColor(kMagenta);
    elZY5->SetLineWidth(1);
    elZY5->Draw();



    TPad * paderZX = (TPad*)can->cd(7);
    TGraph  evdZX = objs[0].f_graph_3D[2];//new TGraph(reco_shower_sp_x.size(), &(reco_shower_sp_z)[0], &(reco_shower_sp_x)[0]);
    TGraph * evdRecoStartPointZX = new TGraph(1,&(reco_vertex_3D)[2], &(reco_vertex_3D)[0]);
    TGraph * evdTrueStartPointZX = new TGraph(1,&(true_vertex)[2], &(true_vertex)[0]);
    evdZX.SetMarkerStyle(20);
    evdZX.SetMarkerSize(1);
    evdZX.SetMarkerColor(objs[0].f_col);
    evdZX.SetTitle("");
    evdZX.GetYaxis()->SetTitle(" X [cm]");
    evdZX.GetXaxis()->SetTitle(" Z [cm]");
    evdZX.GetYaxis()->SetRangeUser((1.0-buffer)*X_min,(1.0+buffer)*X_max);
    evdZX.GetXaxis()->SetLimits((1.0-buffer)*Z_min, (1.0+buffer)*Z_max);
    evdZX.Draw("ap");

    for(int i=1; i<objs.size(); i++){
        objs[i].f_graph_3D[2].SetMarkerStyle(20);
        objs[i].f_graph_3D[2].SetMarkerSize(1);
        objs[i].f_graph_3D[2].SetMarkerColor(objs[i].f_col);
        objs[i].f_graph_3D[2].Draw("p same");
    }

    evdRecoStartPointZX->Draw("same p");
    evdRecoStartPointZX->SetMarkerStyle(29);
    evdRecoStartPointZX->SetMarkerSize(3);
    evdRecoStartPointZX->SetMarkerColor(kRed);
    evdTrueStartPointZX->Draw("same p");
    evdTrueStartPointZX->SetMarkerStyle(29);
    evdTrueStartPointZX->SetMarkerSize(2);
    evdTrueStartPointZX->SetMarkerColor(kCyan);

    TEllipse *elZX3 = new TEllipse( reco_vertex_3D[2],reco_vertex_3D[0],3,3);
    elZX3->SetFillStyle(0);
    elZX3->SetLineColor(kRed-6);
    elZX3->SetLineWidth(1);
    elZX3->Draw();
    TEllipse *elZX5 = new TEllipse( reco_vertex_3D[2],reco_vertex_3D[0],radius,radius);
    elZX5->SetFillStyle(0);
    elZX5->SetLineColor(kMagenta);
    elZX5->SetLineWidth(1);
    elZX5->Draw();

    can->cd(8);    
    TLegend *l3d = new TLegend(0.11,0.11,0.89,0.89);
    l3d->SetFillStyle(0);
    l3d->Draw();
    l3d->AddEntry(&evdXY,"Reco 3D Spacepoints","p");
    l3d->AddEntry(evdRecoStartPointXY,"Reco 3D Start Point","p");
    l3d->AddEntry(evdTrueStartPointXY,"True 3D Start Point","p");
    l3d->AddEntry(elZX3,"3cm Radius","l");
    l3d->AddEntry(elZX5,(std::to_string(radius)+"cm Radius").c_str(),"l");


    //Bottom 3 pads, is also where we calculate the actual numbers

    //Look at just the first Ncm now


    std::vector<double> zp = {0};
    TGraph zero(1,&zp[0],&zp[0]);
    zero.SetMarkerStyle(29);
    zero.SetMarkerSize(3);
    zero.SetMarkerColor(kBlack);

    if(out_graphs2D.size()!=0){
        for(int i=0; i<3; i++){
            can->cd(9+i);

            std::cout<<"Int Line "<<i<<std::endl;
            out_graphs2D[i].Draw("ap");
            out_graphs2D[i].GetYaxis()->SetRangeUser(-radius,radius);
            out_graphs2D[i].GetXaxis()->SetLimits(-radius,radius);
            out_graphs2D[i].SetTitle("");
            out_graphs2D[i].SetMarkerSize(2);
            out_graphs2D[i].SetMarkerStyle(20);
            out_graphs2D[i].SetMarkerColor(kRed-6);
            out_graphs2D[i+3].Draw("p");
            out_graphs2D[i+3].SetMarkerSize(2);
            out_graphs2D[i+3].SetMarkerStyle(20);
            out_graphs2D[i+3].SetMarkerColor(kGreen-6);

            getLine2D(&left_fit[0],kRed-6,i,radius);
            getLine2D(&right_fit[0],kGreen-6,i,radius);
            zero.Draw("same p");
        }
    }

    //**************************** INFO ***************************/
    TPad *p_top_info = (TPad*)can->cd(4);
    p_top_info->cd();

    std::cout<<"Printing INFO into canvas 4"<<std::endl;
    TLegend l_top(0.1,0.0,0.9,1.0);
    l_top.SetTextSize(0.05);

    TLine *l = new TLine(0,0,0,0);
    l->SetLineColor(kWhite);
    l_top.AddEntry(l,("Reco  #theta e^{+}e^{-} : "+std::to_string(reco_ang)).c_str(),"l");
    for(int p=0; p<tags.size(); p++){

        l_top.AddEntry(l,(tags[p]+" : "+vals[p]).c_str(),"l");
        /*
           if(vec_graphs[p][0].GetN()>0){
           l_top.AddEntry(&vec_graphs[p][0],vec_pfp_legend[p].c_str(),"f");
           }else if(vec_graphs[p][1].GetN()>0){
           l_top.AddEntry(&vec_graphs[p][1],vec_pfp_legend[p].c_str(),"f");
           }else if(vec_graphs[p][2].GetN()>0){
           l_top.AddEntry(&vec_graphs[p][2],vec_pfp_legend[p].c_str(),"f");
           }

*/
    }

    l_top.SetHeader((print_name+". 2D recob::Hits").c_str(),"C");
    l_top.SetLineWidth(0);
    l_top.SetLineColor(kWhite);
    l_top.Draw("same");

    std::cout<<"And save PDF"<<std::endl;
    can->Update();
    can->SaveAs((print_name+".pdf").c_str(),"pdf");
    can->Close(); gSystem->ProcessEvents();
    return;
	}//End of plotter

};

