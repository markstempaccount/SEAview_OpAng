#ifndef RECOOPANG1_H
#define RECOOPANG1_H

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
#include <cmath>
#include <iostream>
#include <ostream>
#include <fstream>
#include <map>
#include "TRandom3.h"
#include "TLine.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"


//Some useful functions. Distance from point to plane in 3D. Plane parameterized by ax+by+cz+d =0
 double pointPlaneDistance(double x,double y,double z, double a, double b, double c);

 
// define the parametric line equation
 void line(double R, const double *p, double &x, double &y, double &z) ;

 
 // function Object to be minimized. The distance from a line to a scatter of points squared
 struct SumDistance2 {
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
 
    bool first;
    SumDistance2(TGraph2D *g) : fGraph(g) {first=true;}
 
    // calculate distance line-point
    double distance2(double x,double y,double z, const double *p) {
       // distance line point is D= | (xp-x0) cross  ux |
       // where ux is direction of line and x0 is a point in the line (like t = 0)
        ROOT::Math::XYZVector xp(x,y,z);
        ROOT::Math::XYZVector x0 (0.0, 0.0, 0.);
        ROOT::Math::XYZVector x1(5.0*cos(p[0])*sin(p[1]), 5.0*sin(p[0])*sin(p[1]), 5.0*cos(p[1]) );
        ROOT::Math::XYZVector u = (x1-x0).Unit();
       double d2 = ((xp-x0).Cross(u)).Mag2();
       return d2;
    }
 
    // implementation of the function to be minimized
    double operator() (const double *par) {
       assert(fGraph != 0);
       double * x = fGraph->GetX();
       double * y = fGraph->GetY();
       double * z = fGraph->GetZ();
       int npoints = fGraph->GetN();
       double sum = 0;
       for (int i  = 0; i < npoints; ++i) {
          double d = distance2(x[i],y[i],z[i],par);
          sum += d;
       }
       if (first) {
          std::cout << "Total Initial distance square = " << sum << std::endl;
       }
       first = false;
       return sum;
    }
 
 };
 
struct SumPlaneDistance2 {
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
    bool first;
    SumPlaneDistance2(TGraph2D *g) : fGraph(g) {first=true;}
 
    // calculate distance line-point
    double distance2(double x,double y,double z, const double *p) {
       return (p[0]*x+p[1]*y+p[2]*z)/sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    }
 
    // implementation of the function to be minimized
    double operator() (const double *par) {
       assert(fGraph != 0);
       double * x = fGraph->GetX();
       double * y = fGraph->GetY();
       double * z = fGraph->GetZ();
       int npoints = fGraph->GetN();
       double sum = 0;
       for (int i  = 0; i < npoints; ++i) {
          double d = distance2(x[i],y[i],z[i],par);
          sum += d*d;
       }
       if (first) {
          std::cout << "Total Initial distance square = " << sum << std::endl;
       }
       first = false;
       return sum;
    }
 
 };
 
//Given input parameters plot line in various 2Ds
int getLine2D(double * parFit, int col,int which,double radius);

int getLine(const double * parFit, int col);

//Messy at the moment. Will need to tidy up
double recoOpAng1(std::vector<double> start_point,   std::vector<double> shower_point_x,     std::vector<double> shower_point_y,    std::vector<double> shower_point_z, std::vector<TGraph> &out_graphs2D, std::vector<double> &left_par, std::vector<double>&right_par);
#endif
