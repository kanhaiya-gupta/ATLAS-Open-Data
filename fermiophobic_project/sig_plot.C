#include <TH1.h>
#include <TF2.h>
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
#include "TGraph.h"

void sig_plot(){

  gStyle->SetOptStat(0);
  //  (new TCanvas())->DrawFrame(0.1, 0.1, 2.1, 2.1);
   TFile *Data = TFile::Open("data.root");


   TH2D* data  = (TH2D*) Data->Get("lepjet2d");

   TEllipse *el = new TEllipse(1.,1.,.037,.576);
  
  
   // TF2 *f2 = new TF2("f2","((x-1)/0.0001)^2 + ((y-1)/0.0001)^2 -1.",0.,2.,0.,2.);  // PDF to be estimated
   // TLine *line = new TLine(-3,0,8,0);
   data->SetStats(0);
   data->GetXaxis()->SetTitle("M_{ll}/92 (GeV/c^{2})");
   data->GetYaxis()->SetTitle("M_{bb}/125 (GeV/c^{2})");
   data->Draw("colz");
   el->SetFillStyle(0);
   el->SetLineWidth(5);
   el->SetLineColor(kRed);
   el->Draw("same");

   //  TH1D* data1  = (TH1D*) Data->Get("dibjet");
   // data1->Draw();
  
   TLegend *legend=new TLegend(0.6,0.65,0.78,0.80);    // x1, y1, x2, y2 coordinates of the Legend 
   // legend->SetTextFont(72);
   legend->SetTextSize(0.04);
   legend->SetLineColor(0);
   //  legend->AddEntry(data,"CL","l");
   legend->AddEntry(el,"el","Signal Region ");
   //  legend->Draw();

  

  
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

   
