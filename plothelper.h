#ifndef PLOTHELPER_H
#define PLOTHELPER_H

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
#include <iomanip>      // std::setprecision

#include <iostream>
#include <fstream>
#include <map>
#include "TRandom3.h"
#include "TLine.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"

#include <string>


template <typename T>
std::string to_string_prec(const T a_value, const int n)
{
    std::ostringstream out;
    out <<std::fixed<< std::setprecision(n) << a_value;
    return out.str();
}
// -----------  Functions to get a TH1 given a TTree, variable, cuts, uniq name, and binning (in usual (Nbins, min,max) form) -------- 
TH1 * getTH1(TTree * tin, std::string var, std::string cuts, std::string nam, std::string binning);

// -----------  Functions to get a TH1 given a TTree, variable, cuts, uniq name, and binning (in a form of a vector of bin boundaries) -------- 
TH1D *getTH1(TTree *tin, std::string var, std::string cuts, std::string nam, std::vector<double>  bins);
TTree* loadgLEE(std::string filename, std::string int_dir);


double dist2D(double w1, double t1, double w2, double t2);


#endif
