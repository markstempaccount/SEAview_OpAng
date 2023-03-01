#include "SetRealAspectRatio2.h"

bool SetRealAspectRatio2(TCanvas *ch, const Int_t axis)
{
   ch -> Update();
 
   //Get how many pixels are occupied by the canvas
   Int_t npx = ch -> GetWw();
   Int_t npy = ch -> GetWh();
 
   //Get x-y coordinates at the edges of the canvas (extrapolating outside the axes, NOT at the edges of the histogram)
   Double_t x1 = ch -> GetX1();
   Double_t y1 = ch -> GetY1();
   Double_t x2 = ch -> GetX2();
   Double_t y2 = ch -> GetY2();
 
   //Get the length of extrapolated x and y axes
   Double_t xlength2 = x2 - x1;
   Double_t ylength2 = y2 - y1;
   Double_t ratio2   = xlength2/ylength2;
 
   //Now get the number of pixels including the canvas borders
   Int_t bnpx = ch -> GetWindowWidth();
   Int_t bnpy = ch -> GetWindowHeight();
 
   if (axis==1) {
      ch -> SetCanvasSize(TMath::Nint(npy*ratio2), npy);
      ch -> SetWindowSize((bnpx-npx)+TMath::Nint(npy*ratio2), bnpy);
   } else if (axis==2) {
      ch -> SetCanvasSize(npx, TMath::Nint(npx/ratio2));
      ch -> SetWindowSize(bnpx, (bnpy-npy)+TMath::Nint(npx/ratio2));
   } else {
      Error("SetRealAspectRatio", "axis value %d is neither 1 (resize along x-axis) nor 2 (resize along y-axis).",axis);
      return false;
   }
 
   //Check now that resizing has worked
 
   ch -> Update();
 
   //Get how many pixels are occupied by the canvas
   npx = ch -> GetWw();
   npy = ch -> GetWh();
 
   //Get x-y coordinates at the edges of the canvas (extrapolating outside the axes,
   //NOT at the edges of the histogram)
   x1 = ch -> GetX1();
   y1 = ch -> GetY1();
   x2 = ch -> GetX2();
   y2 = ch -> GetY2();
 
   //Get the length of extrapolated x and y axes
   xlength2 = x2 - x1;
   ylength2 = y2 - y1;
   ratio2 = xlength2/ylength2;
 
   //Check accuracy +/-1 pixel due to rounding
   if (abs(TMath::Nint(npy*ratio2) - npx)<2) {
      return true;
   } else {
      Error("SetRealAspectRatio", "Resizing failed.");
      return false;
   }
}

