#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
using std::vector;
#include <iostream>
#include <fstream>
#include "TChain.h"
#include "TMath.h"
#include "TObject.h"
#include <TVirtualFitter.h>


void new_fit(){
  TFile *fileData = TFile::Open("data.root");  //sc_data.root
  TFile *fileMC1 = TFile::Open("WZ.root");   // sc2_WZ.root
  TFile *fileMC2 = TFile::Open("ZZ.root"); // sc2_ZZ.root

   TH1F* data = (TH1F*) fileData->Get("M_T(W)");    // data histogram
   TH1F* mc1 = (TH1F*) fileMC1->Get("M_T(W)");     // first weighted  MC histogram
   TH1F* mc2 = (TH1F*) fileMC2->Get("M_T(W)");     // second weighted MC histogram

   // TH1F* mc3 = (TH1F*) fileMC1->Get("M_T_weit_mc(W)");     // M_T_weit(W)->Divide(M_T_mc(W)), M_T_mc(W) unweighted for mc1
   // TH1F* mc4 = (TH1F*) fileMC2->Get("M_T_weit_mc(W)");     // similarly for zz

   

  
 
   TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
   mc->Add(mc1);
   mc->Add(mc2);
   TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  
   
   fit->Constrain(0,0.0,1.0);               // constrain fraction 0 to be between 0 and 1
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->SetRangeX(1,40);                    // use only from the 7th bin to 22th that has non zero entry.
   //  fit->SetWeight(0,mc3);
   // fit->SetWeight(1,mc4);
   
   Int_t status = fit->Fit();               // perform the fit
   std::cout << "fit status: " << status << std::endl;
   if (status == 0) {                       // check on fit status
     TH1F* result = (TH1F*) fit->GetPlot();
     data->Draw("Ep");
     result->Draw("same");
   }

   Double_t  value0, error0, value1, error1;
  fit->GetResult(0, value0, error0);
  fit->GetResult(1, value1, error1);

  //  TVirtualFitter* vFit = fit->GetFitter();
}
