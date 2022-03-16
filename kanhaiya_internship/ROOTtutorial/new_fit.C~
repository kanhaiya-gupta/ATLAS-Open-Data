#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
using std::vector;
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"

void new_fit(){
   TFile *fileData = TFile::Open("DataFull.root");
   TFile *fileMC1 = TFile::Open("WZ.root");
   TFile *fileMC2 = TFile::Open("ZZ.root");

   TH1D* data = (TH1D*) fileData->Get("M_T(W)");    // data histogram
   TH1D* mc1 = (TH1D*) fileMC1->Get("M_T(W)");     // first MC histogram
   TH1D* mc2 = (TH1D*) fileMC2->Get("M_T(W)");     // second MC histogram

   // Scaling;

   mc1->Scale(0.79);
   mc2->Scale(0.21);
 
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   mc->Add(mc1);
   mc->Add(mc2);
   TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(1,100);                    // use only the first 50 bins in the fit
   Int_t status = fit->Fit();               // perform the fit
   std::cout << "fit status: " << status << std::endl;
   if (status == 0) {                       // check on fit status
     TH1D* result = (TH1D*) fit->GetPlot();
     data->Draw("Ep");
     result->Draw("same");
   }
}
