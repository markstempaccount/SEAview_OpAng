#ifndef PLOTHELPER_H
#define PLOTHELPER_H


template <typename T>
std::string to_string_prec(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out <<std::fixed<< std::setprecision(n) << a_value;
    return out.str();
}

// -----------  Functions to get a TH1 given a TTree, variable, cuts, uniq name, and binning (in usual (Nbins, min,max) form) -------- 
TH1 * getTH1(TTree * tin, std::string var, std::string cuts, std::string nam, std::string binning){

    tin->Draw((var+">>"+nam+ binning).c_str() , ("("+cuts+ ")").c_str(),"goff");
    TH1* th1 = (TH1*)gDirectory->Get(nam.c_str()) ;
    th1->SetLineWidth(1);
    th1->SetStats(0);
    th1->SetDirectory(0);

    return th1;
}


// -----------  Functions to get a TH1 given a TTree, variable, cuts, uniq name, and binning (in a form of a vector of bin boundaries) -------- 
TH1D *getTH1(TTree *tin, std::string var, std::string cuts, std::string nam, std::vector<double>  bins){

    TH1D* th1 = new TH1D(nam.c_str(),nam.c_str(), bins.size()-1, &(bins)[0]);
    tin->Draw((var+">>+"+nam).c_str() , ("("+cuts+")").c_str(),"goff");
    th1->SetLineWidth(1);
    th1->SetStats(0);
    th1->SetDirectory(0);

    return th1;

}

TTree* loadgLEE(std::string filename, std::string int_dir){

    TFile*f = new TFile((filename).c_str(),"read");
    TTree*v = (TTree*)f->Get((int_dir+"/vertex_tree").c_str());
    TTree*s = (TTree*)f->Get((int_dir+"/simple_tree").c_str());
    TTree*e = (TTree*)f->Get((int_dir+"/eventweight_tree").c_str());
    v->AddFriend(s);
    v->AddFriend(e);

    return v;
}

double dist2D(double w1, double t1, double w2, double t2){
    double wire_con=0.3;
    double tick_con=1.0/25.0;
    return  sqrt(pow(w1*wire_con-w2*wire_con,2)+pow(t1*tick_con-t2*tick_con,2));
}


#endif
