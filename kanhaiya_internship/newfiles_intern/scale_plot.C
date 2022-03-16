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
#include <cstdlib> 
#include <cstdio> 
#include <string> 
#include <iostream> 
#include <fstream>
#include "TSystem.h"
#include "myroot.h"
#include "TGraph.h"

void scale_plot(){

   TFile *Data = TFile::Open("data.root");
   TFile *fileMC1 = TFile::Open("WZ.root");
   TFile *fileMC2 = TFile::Open("ZZ.root");
   TFile *fileMC3 = TFile::Open("BSM1.root");

   TH1D* data  = (TH1D*) Data->Get("mass_trimu_sig");
   TH1D* mc1 = (TH1D*) fileMC1->Get("mass_trimu_sig");
   TH1D* mc2 = (TH1D*) fileMC2->Get("mass_trimu_sig");
   TH1D* bsm = (TH1D*) fileMC3->Get("mass_trimu_sig");

   float norm = 1.;
   float scale1 = norm/data->Integral();
   data->Scale(scale1,"width");

   float scale2 = norm/mc1->Integral();
   mc1->Scale(scale2,"width");

   float scale3 = norm/mc2->Integral();
   mc2->Scale(scale3,"width");

   float scale4 = norm/bsm->Integral();
   bsm->Scale(scale4,"width");
   
   float n_obs  = data->Integral();
   float n_bkg = mc1->Integral() +  mc2->Integral();
   float n_sig = bsm->Integral();

   std::cout<<"Number of observed events: "<<n_obs <<std::endl;
   std::cout<<"Number of background events: "<<n_bkg <<std::endl;
   std::cout<<"Number of expected signal: "<<n_sig <<std::endl;

   float lumi_int = 1.0;
   float cl = 0.95;

   long int nsamples = 100;
   vector<float> CL_int;
   vector<float> sig_int;

   for (long int i=0; i<nsamples; i++)
  {
    float sig = i*1.0;
   float mu = ( sig + n_bkg)*lumi_int;
   float beta = ROOT::Math::poisson_cdf(n_obs,mu)/ROOT::Math::poisson_cdf(n_obs,n_bkg);
   float CL = 1 - beta;

   CL_int.push_back(CL);
   sig_int.push_back(sig);
  
  }

   float s_up = 0.5*TMath::ChisquareQuantile(cl,2*(n_obs+1))/lumi_int;
   TGraph* gr = new TGraph(nsamples,&sig_int[0],&CL_int[0]);
   TF1 *f1 = new TF1("f1", "0.95", 0,100);  // PDF to be estimated
   // TLine *line = new TLine(-3,0,8,0);
   gr->SetLineWidth(2);
   gr->Draw("");
   f1->Draw("same");
   gr->SetTitle("Confidence Interval");
   gr->GetYaxis()->SetTitle("CL");
   gr->GetXaxis()->SetTitle("#epsilon #sigma");

   TLegend *legend=new TLegend(0.6,0.65,0.78,0.80);    // x1, y1, x2, y2 coordinates of the Legend 
   // legend->SetTextFont(72);
   legend->SetTextSize(0.04);
   legend->SetLineColor(0);
   legend->AddEntry(gr,"CL, n_{obs} = 0","l");
   legend->AddEntry(f1,"CL=95%","l");
   legend->Draw();

   //  std::cout<<"Confidence Interval at 90% : "<<CL <<std::endl;   
   std::cout<<"critical value: "<<s_up <<std::endl;   

  

  
   //  Math::poisson_cdf
   
  /*
   auto c = new TCanvas("c", "c", 600, 600);
   c->Divide(1,2);
   TFile *fileMC1 = TFile::Open("sc_WZ.root");
   TFile *fileMC2 = TFile::Open("sc_ZZ.root");

  
   TH1D* mc1 = (TH1D*) fileMC1->Get("M_T(W)");     // first MC histogram
   TH1D* mc2 = (TH1D*) fileMC1->Get("M_T_s(W)");     // second MC histogram
   TH1D* mc3 = (TH1D*) fileMC1->Get("M_T(WZ)");
   TH1D* mc4 = (TH1D*) fileMC1->Get("M_T_s(WZ)");

   TH1D* mc5 = new TH1D(*mc2);
   TH1D* mc6 = new TH1D(*mc4);

   mc5->Divide(mc1);
   mc6->Divide(mc3);

   c->cd(1);
   mc5->Draw("Hist");

   c->cd(2);
   mc6->Draw("Hist");
   
   
   c->cd(1);
   mc1->SetMarkerStyle(34);
   mc1->Draw();
   mc2->SetMarkerStyle(4);
   mc2->Draw("same C");

   c->cd(2);
   mc3->SetMarkerStyle(34);
   mc3->Draw();
   mc4->SetMarkerStyle(4);
   mc4->Draw("same C");

   float WZ_W = mc1->Integral();
   float ZZ_W = mc2->Integral();
   float WZ_NW = mc3->Integral();
   float ZZ_NW = mc4->Integral();

   float weight_WZ = WZ_W/WZ_NW;
   float weight_ZZ = ZZ_W/ZZ_NW;

   std::cout<<"Weight for WZ: "<<weight_WZ<<std::endl;
   std::cout<<"Weight for ZZ: "<<weight_ZZ<<std::endl;
   */
}

   
