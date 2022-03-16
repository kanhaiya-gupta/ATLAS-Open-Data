#include "myroot.h"

double fitf(Double_t *x, Double_t *par)
{
   double fitval = par[0]*TMath::Exp(par[1]*x[0]);
   return fitval;
}

void fitExample()
{
   gStyle->SetOptFit(11111);
   double xmin=0.04, xmax=0.4; // fit range
   int npara = 2;  // number of parameters to be fitted
   TString hname="example";
   
   // obtaining the histogram
   TH1D* hist = (TH1D*) gROOT->FindObject(hname);
   if (!hist){
      cout << "opening file histo1.root" << endl;
      TFile *file = TFile::Open("histo1.root");
      hist = (TH1D*) file->Get(hname);
   }
   hist->SetMarkerStyle(20);
   hist->SetMarkerSize(0.8);
   hist->SetMarkerColor(kBlue);
   
   TCanvas *c1 = new TCanvas("c1","the fit canvas",500,400);
// Creates a Root function based on function fitf above
   TF1 *func = new TF1("fitf",fitf,xmin,xmax,npara);
// Sets initial values and parameter names
   func->SetParameters(100.,-1.); // initial parameters before minimization
   func->SetParNames("Constant","tau");
// Fit histogram in range defined by function
   hist->Fit(func,"r,e"); // r: fit in function range, e: use MINOS for error calc.
   hist->Draw("e,same");

// Get integral of function between fit limits
   printf("Integral of function = %g\n",func->Integral(xmin, xmax));
}
