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

   gStyle->SetOptStat(0);
   TFile *fileData = TFile::Open("data.root");
   TFile *fileMC1 = TFile::Open("Zjet.root");
   TFile *fileMC2 = TFile::Open("ZZ.root");
   TFile *fileMC3 = TFile::Open("ttbar.root");
   TFile *fileMC4 = TFile::Open("WZ.root");
   

   TH1F* data = (TH1F*) fileData->Get("lepjet");    // data histogram
   TH1F* mc1 = (TH1F*) fileMC1->Get("lepjet");    // first weighted  MC histogram
   TH1F* mc2 = (TH1F*) fileMC2->Get("lepjet");     // second weighted MC histogram
   TH1F* mc3 = (TH1F*) fileMC3->Get("lepjet");     // third weighted MC histogram
   TH1F* mc4 = (TH1F*) fileMC4->Get("lepjet");     // fourth MC histogram
   


  
 
   TObjArray *mc = new TObjArray(4);        // MC histograms are put in this array
   mc->Add(mc1);
   mc->Add(mc2);
   mc->Add(mc3);
   mc->Add(mc4);
   TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
  
   
   fit->Constrain(0,0.0,1.0);               // constrain fraction 0 to be between 0 and 1
   fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->Constrain(2,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
   fit->Constrain(3,0.0,1.0);
   //  fit->SetRangeX(3,15);                    // use only from the 7th bin to 22th that has non zero entry.
   
  
   
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
  //  fit->GetResult(2, value2, error2);

  //  TVirtualFitter* vFit = fit->GetFitter();
}
