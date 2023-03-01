#ifndef SETREALASPECTRATIO2_H
#define SETREALASPECTRATIO2_H

#include <TCanvas.h>
#include <TMath.h>
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

#include <string>


bool SetRealAspectRatio2(TCanvas *ch, const Int_t axis);

#endif
