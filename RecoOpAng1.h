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
 

//Some useful functions. Distance from point to plane in 3D. Plane parameterized by ax+by+cz+d =0
 double pointPlaneDistance(double x,double y,double z, double a, double b, double c) {
       return (a*x+b*y+c*z)/sqrt(a*a+b*b+c*c);
    }

 
// define the parametric line equation
 void line(double R, const double *p, double &x, double &y, double &z) {
    // a parametric line in spherical polar coordinates is defined by 2 parameters, theta and phi and moving along the radius R
    x = R*cos(p[0])*sin(p[1]);   //moving to cartesian
    y = R*sin(p[0])*sin(p[1]);
    z = R*cos(p[1]);
 }
 
 bool first = true;
 
 // function Object to be minimized. The distance from a line to a scatter of points squared
 struct SumDistance2 {
    // the TGraph is a data member of the object
    TGraph2D *fGraph;
 
    SumDistance2(TGraph2D *g) : fGraph(g) {}
 
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
 
    SumPlaneDistance2(TGraph2D *g) : fGraph(g) {}
 
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
int getLine2D(double * parFit, int col,int which,double radius){
    std::vector<double> parFit2 = {parFit[0],parFit[1]};
    parFit2[1]=parFit2[1]-TMath::Pi();

    int n = 1000;
    double t0 = 0.0;
    double dt = radius;
    TPolyLine *l = new TPolyLine(n);
    TPolyLine *l2 = new TPolyLine(n);
    for (int i = -n; i <n;++i) {
        double t = t0+ dt*i/n;
        double x,y,z;
        double xval,yval;
        
        line(t,parFit,x,y,z);
        switch(which){
            case(0): 
                xval = x; yval = y;
                break;
            case(1):
                xval = z; yval = y;
                break;
            case(2):
                xval =z;yval =x;
                break;
            default:
                xval = z;yval=y;
                break;
        }

        l->SetPoint(i,xval,yval);
        
        line(t,&(parFit2)[0],x,y,z);
        switch(which){
            case(0): 
                xval = x; yval = y;
                break;
            case(1):
                xval = z; yval = y;
                break;
            case(2):
                xval =z;yval =x;
                break;
            default:
                xval = z;yval=y;
                break;
        }

        l2->SetPoint(i,xval,yval);
        
    }
    l->SetLineColor(col);
    l->SetLineWidth(1);
    l->Draw("same");
    l2->SetLineColor(col);
    l2->SetLineWidth(1);
    l2->Draw("same");
 return 0;
}
int getLine(const double * parFit, int col){

    std::vector<double> parFit2 = {parFit[0],parFit[1]};
    parFit2[1]=parFit2[1]-TMath::Pi();

    int n = 1000;
    double t0 = 0.0;
    double dt = 5.0;
    TPolyLine3D *l = new TPolyLine3D(n);
    TPolyLine3D *l2 = new TPolyLine3D(n);
    for (int i = -n; i <n;++i) {
        double t = t0+ dt*i/n;
        double x,y,z;
        line(t,parFit,x,y,z);
        l->SetPoint(i,x,y,z);
        line(t,&(parFit2)[0],x,y,z);
        l2->SetPoint(i,x,y,z);
    }
    l->SetLineColor(col);
    l->SetLineWidth(2);
    l->Draw("same");
    l2->SetLineColor(col);
    l2->SetLineWidth(2);
    l2->Draw("same");
 return 0;
}

//Messy at the moment. Will need to tidy up
double recoOpAng1(std::vector<double> start_point,   std::vector<double> shower_point_x,     std::vector<double> shower_point_y,    std::vector<double> shower_point_z, std::vector<TGraph> &out_graphs2D, std::vector<double> &left_par, std::vector<double>&right_par){

    int N_pts = shower_point_x.size();
    
    std::vector<double> null_point = {0,0,0};
    gStyle->SetOptStat(0);
    gStyle->SetOptFit();

    for(auto &x: shower_point_x) x-=start_point[0];
    for(auto &y: shower_point_y) y-=start_point[1];
    for(auto &z: shower_point_z) z-=start_point[2];


    TGraph2D *g3D = new TGraph2D((std::to_string(start_point[0])+std::to_string(shower_point_x[0]*shower_point_y[0])).c_str(),"",N_pts, &(shower_point_x)[0], &(shower_point_y)[0], &(shower_point_z)[0]);
    g3D->SetMarkerSize(2);
    g3D->SetMarkerStyle(20);
    //g3D->Fit("pol1");
    //g3D->Draw("ap"); 

    ROOT::Fit::Fitter  fitter;
    
    // make the functor objet
    SumDistance2 sdist(g3D);
    ROOT::Math::Functor fcn(sdist,2);

    //fitter.GetFCN()->SetParLimit(0,0,2.0*TMath::Pi());
    //fitter.GetFCN()->SetParLimit(1,0.0,TMath::Pi());

    // set the function and the initial parameter values
    double pStart[2] = {0,0};

    fitter.SetFCN(fcn,pStart);
    // set step sizes different than default ones (0.3 times parameter values)
//    for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.1);

//    g3D->Draw("p0");

//Plot the origin
    TGraph2D * origin = new TGraph2D( (std::to_string(start_point[0])+std::to_string(shower_point_x[0]*shower_point_y[0])+"RecpSP").c_str(),"",1,&(null_point)[0], &(null_point)[1], &(null_point)[2]);
//    evdRecoStartPoint->Draw("same p");
    origin->SetMarkerStyle(29);
    origin->SetMarkerSize(5);
    origin->SetMarkerColor(kMagenta);


    bool ok = fitter.FitFCN();
    if (!ok) {
        Error("line3Dfit","Line3D Fit failed");
        return -9995;
    }

    const ROOT::Fit::FitResult & result = fitter.Result();

    std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
    result.Print(std::cout);


    // get fit parameters
    const double * parFit = result.GetParams();
   

      /*
    std::vector<double> parFit2 = {parFit[0],parFit[1]};
    parFit2[1]=parFit2[1]-TMath::Pi();

    // draw the fitted line
    int n = 1000;
    double t0 = 0.0;
    double dt = 5.0;
    TPolyLine3D *l = new TPolyLine3D(n);
    TPolyLine3D *l2 = new TPolyLine3D(n);
    for (int i = -n; i <n;++i) {
        double t = t0+ dt*i/n;
        double x,y,z;
        line(t,parFit,x,y,z);
        l->SetPoint(i,x,y,z);
        line(t,&(parFit2)[0],x,y,z);
        l2->SetPoint(i,x,y,z);
    }
    l->SetLineColor(kBlack);
    l->SetLineWidth(2);
    l->Draw("same");
    l2->SetLineColor(kBlack);
    l2->SetLineWidth(2);
    l2->Draw("same");
        */
    //getLine(parFit,kBlack);


    //Now fit a plane. 

    ROOT::Fit::Fitter  planefitter;
    SumPlaneDistance2 pdist(g3D);
    ROOT::Math::Functor planefcn(pdist,3);
   
    // set the function and the initial parameter values
    double planeStart[3] = {0.1,0.1,0.1};
    planefitter.SetFCN(planefcn,planeStart);

    bool pok = planefitter.FitFCN();
    if (!pok) {
        Error("Plane3Dfit","Plane3D Fit failed");
        return -9996;
    }

    const ROOT::Fit::FitResult & planeresult = planefitter.Result();

    planeresult.Print(std::cout);

    const double * planeparFit = planeresult.GetParams();
    std::cout<<"Plane Normal is : "<<planeparFit[0]<<" , "<<planeparFit[1]<<" , "<<planeparFit[2]<<std::endl;
    std::cout<<"Norm: "<<sqrt( pow(planeparFit[0],2)+pow(planeparFit[1],2)+pow(planeparFit[2],2))<<std::endl;
    std::cout<<"Line is : "<<cos(parFit[0])*sin(parFit[1])<<" , "<<sin(parFit[0])*sin(parFit[1])<<" , "<<cos(parFit[1])<<std::endl   ;

    ROOT::Math::XYZVector xp(planeparFit[0],planeparFit[1],planeparFit[2]);
    ROOT::Math::XYZVector xl (cos(parFit[0])*sin(parFit[1]), sin(parFit[0])*sin(parFit[1]), cos(parFit[1]));
    ROOT::Math::XYZVector res = (xp.Cross(xl));
    std::cout<<"Cross-P "<<res.X()<<" "<<res.Y()<<" "<<res.Z()<<std::endl;



    std::vector<double> left_x;
    std::vector<double> left_y;
    std::vector<double> left_z;

    std::vector<double> right_x;
    std::vector<double> right_y;
    std::vector<double> right_z;

    for(int i=0; i< N_pts; i++){

        double d = pointPlaneDistance(shower_point_x[i],shower_point_y[i], shower_point_z[i], res.X(), res.Y(), res.Z() );
        if(d<=0){
                left_x.push_back(shower_point_x[i]);
                left_y.push_back(shower_point_y[i]);
                left_z.push_back(shower_point_z[i]);
        }else{
                right_x.push_back(shower_point_x[i]);
                right_y.push_back(shower_point_y[i]);
                right_z.push_back(shower_point_z[i]);
        }

    }

    std::cout<<" Left: "<<left_x.size()<<" Right: "<<right_x.size()<<std::endl;
    TGraph2D *left3D = new TGraph2D((std::to_string(start_point[0])+std::to_string(shower_point_x[0]*shower_point_y[0])+"left3D").c_str(),"",left_x.size(), &(left_x)[0], &(left_y)[0], &(left_z)[0]);
    left3D->SetMarkerSize(2);
    left3D->SetMarkerStyle(20);
    left3D->SetMarkerColor(kRed-6);
//    left3D->Draw("p same");

    TGraph2D *right3D = new TGraph2D((std::to_string(start_point[0])+std::to_string(shower_point_x[0]*shower_point_y[0])+"right3D").c_str(),"",right_x.size(), &(right_x)[0], &(right_y)[0], &(right_z)[0]);
    right3D->SetMarkerSize(2);
    right3D->SetMarkerStyle(20);
    right3D->SetMarkerColor(kBlue-6);
//    right3D->Draw("p same");
   
    
    out_graphs2D.emplace_back(left_x.size(), &(left_x)[0], &(left_y)[0]);
    out_graphs2D.emplace_back(left_x.size(), &(left_z)[0], &(left_y)[0]);
    out_graphs2D.emplace_back(left_x.size(), &(left_z)[0], &(left_x)[0]);
    out_graphs2D.emplace_back(right_x.size(), &(right_x)[0], &(right_y)[0]);
    out_graphs2D.emplace_back(right_x.size(), &(right_z)[0], &(right_y)[0]);
    out_graphs2D.emplace_back(right_x.size(), &(right_z)[0], &(right_x)[0]);

 

    ROOT::Fit::Fitter  leftfitter;
    SumDistance2 leftdist(left3D);
    ROOT::Math::Functor leftfcn(leftdist,2);
    double leftStart[2] = {0,0};
    leftfitter.SetFCN(leftfcn,leftStart);
    bool leftok = leftfitter.FitFCN();
    if (!leftok) {
        Error("LEFTfit","LEFT Line3D Fit failed");
        return -9997;
    }
    const ROOT::Fit::FitResult & leftresult = leftfitter.Result();
    const double * leftparFit = leftresult.GetParams();
    //getLine(leftparFit,kRed-7);
    left_par = {leftparFit[0],leftparFit[1]};

    ROOT::Fit::Fitter  rightfitter;
    SumDistance2 rightdist(right3D);
    ROOT::Math::Functor rightfcn(rightdist,2);
    double rightStart[2] = {0,0};
    rightfitter.SetFCN(rightfcn,rightStart);
    bool rightok = rightfitter.FitFCN();
    if (!rightok) {
        Error("RIGHTfit","RIGHT Line3D Fit failed");
        return -9998;
    }
    const ROOT::Fit::FitResult & rightresult = rightfitter.Result();
    const double * rightparFit = rightresult.GetParams();
    //getLine(rightparFit,kBlue-7);
    right_par = {rightparFit[0],rightparFit[1]};

    std::vector<double> leftres  = {cos(leftparFit[0])*sin(leftparFit[1]) , sin(leftparFit[0])*sin(leftparFit[1]) ,cos(leftparFit[1]) } ;
    std::vector<double> rightres  = {cos(rightparFit[0])*sin(rightparFit[1]) , sin(rightparFit[0])*sin(rightparFit[1]) ,cos(rightparFit[1]) } ;


    double left_norm = sqrt(pow(leftres[0],2)+pow(leftres[1],2)+pow(leftres[2],2));
    double right_norm = sqrt(pow(rightres[0],2)+pow(rightres[1],2)+pow(rightres[2],2));
    std::cout<<left_norm<<" "<<right_norm<<std::endl;
   
    double reco_ang = acos(leftres[0]*rightres[0]+ leftres[1]*rightres[1]+leftres[2]*rightres[2])*180.0/TMath::Pi()/(left_norm*right_norm);

    if(reco_ang > 90) reco_ang = 180.0-reco_ang;
    std::cout<<"THE FINAL ANSWER is : "<<reco_ang<<std::endl;

    return reco_ang;

}
#endif
