#ifndef SEAOBJECT_H
#define SEAOBJECT_H
#include "TH1D.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TEllipse.h"
#include "TLegend.h"

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



class SEAobject{

    public:
    int f_type;
    int f_col;
    std::string f_tag;

    std::vector<double> f_sp_x;
    std::vector<double> f_sp_y;
    std::vector<double> f_sp_z;
    std::vector<double> f_weights;
    int f_num_sp;
    std::vector<double> f_max_sp;
    std::vector<double> f_min_sp;
    
    std::vector<std::vector<double>> f_ticks;
    std::vector<std::vector<double>> f_wires;
    std::vector<std::vector<double>> f_energies;

    double f_max_tick;
    double f_min_tick;

    std::vector<double> f_min_wires;
    std::vector<double> f_max_wires;

    std::vector<size_t> f_num_hits;

    std::vector<TGraph> f_graph_2D;
    std::vector<TGraph> f_graph_3D;

    SEAobject(int type, int col, std::vector<int> wires, std::vector<int> planes, std::vector<double> ticks, std::vector<double> energies, std::vector<double> sp_x,  std::vector<double> sp_y, std::vector<double> sp_z, std::vector<double> weights) : f_type(type), f_col(col), f_sp_x(sp_x), f_sp_y(sp_y), f_sp_z(sp_z), f_weights(weights) {


        f_ticks.resize(3);
        f_wires.resize(3);
        f_energies.resize(3);

        f_num_sp = sp_x.size();

        f_tag = std::to_string(type)+"_"+std::to_string(col)+"_"+std::to_string(f_num_sp);

        int n_hits = planes.size();
        for(int i=0; i< n_hits; i++){
            f_ticks[planes.at(i)].push_back(ticks.at(i)); 
            f_wires[planes.at(i)].push_back(wires.at(i)); 
            f_energies[planes.at(i)].push_back(energies.at(i)); 
        }

        f_num_hits = {f_ticks[0].size(),f_ticks[1].size(), f_ticks[2].size()};
        
        f_max_tick = *max_element(ticks.begin(),ticks.end());
        f_min_tick = *min_element(ticks.begin(), ticks.end());



        f_max_wires = { 
            f_wires[0].size()!=0 ? *max_element(f_wires[0].begin(),f_wires[0].end()) : -9999 ,
            f_wires[1].size()!=0 ? *max_element(f_wires[1].begin(),f_wires[1].end()) : -9999, 
            f_wires[2].size()!=0 ? *max_element(f_wires[2].begin(),f_wires[2].end()) : -9999 }; 
        f_min_wires = { 
            f_wires[0].size()!=0 ? *min_element(f_wires[0].begin(),f_wires[0].end()) : 9999 ,
            f_wires[1].size()!=0 ? *min_element(f_wires[1].begin(),f_wires[1].end()) : 9999, 
            f_wires[2].size()!=0 ? *min_element(f_wires[2].begin(),f_wires[2].end()) : 9999 };

        
        f_graph_2D.emplace_back((int)f_num_hits[0],&(f_wires[0])[0],&(f_ticks[0])[0]);
        f_graph_2D.emplace_back((int)f_num_hits[1],&(f_wires[1])[0],&(f_ticks[1])[0]);
        f_graph_2D.emplace_back((int)f_num_hits[2],&(f_wires[2])[0],&(f_ticks[2])[0]);

        f_graph_3D.emplace_back((int)f_num_sp,&(f_sp_x)[0],&(f_sp_y)[0]);
        f_graph_3D.emplace_back((int)f_num_sp,&(f_sp_z)[0],&(f_sp_y)[0]);
        f_graph_3D.emplace_back((int)f_num_sp,&(f_sp_z)[0],&(f_sp_x)[0]);

        f_max_sp = { *max_element(f_sp_x.begin(),f_sp_x.end()),  *max_element(f_sp_y.begin(),f_sp_y.end()), *max_element(f_sp_z.begin(),f_sp_z.end())};  
        f_min_sp = { *min_element(f_sp_x.begin(),f_sp_x.end()),  *min_element(f_sp_y.begin(),f_sp_y.end()), *min_element(f_sp_z.begin(),f_sp_z.end())};  


    };

    int print(){
        std::cout<<"SEAobject of type "<<f_type<<" has "<<f_num_sp<<" spacepoints and "<<f_num_hits[0]<<" , "<<f_num_hits[1]<<" , "<<f_num_hits[2]<<" hits on each plane. "<<std::endl;
        std::cout<<"SEAobject has min/max ticks of  "<<f_min_tick<<" / "<<f_max_tick<<std::endl;
        std::cout<<"SEAobject has min/max wires of  "<<f_min_wires[0]<<" / "<<f_max_wires[0]<<" on plane0"<<std::endl;
        std::cout<<"SEAobject has min/max wires of  "<<f_min_wires[1]<<" / "<<f_max_wires[1]<<" on plane1"<<std::endl;
        std::cout<<"SEAobject has min/max wires of  "<<f_min_wires[2]<<" / "<<f_max_wires[2]<<" on plane2"<<std::endl;
        return 0;
    };
};



#endif
