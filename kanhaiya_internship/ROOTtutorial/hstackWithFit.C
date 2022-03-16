#include "myroot.h"
#include "FitGausPlusLinear.C"

TCanvas *hstackWithFit() {
// Example of stacked histograms: class THStack
// based on Rene Bruns hstack.C 
   gStyle->SetOptStat(1000000011);
   gStyle->SetOptFit(1111);
   THStack *hs = new THStack("hs","Stacked 1D histograms"); 
   //create four 1-d histograms
   TH1D *h1st = new TH1D("h1st","test hstack",100,-4,4);
   h1st->FillRandom("gaus",20000);
   h1st->SetFillColor(kRed);
   h1st->SetMarkerStyle(21);
   h1st->SetMarkerColor(kRed);
   hs->Add(h1st);
   TH1D *h2st = new TH1D("h2st","test hstack",100,-4,4);
   h2st->FillRandom("gaus",15000);
   h2st->SetFillColor(kBlue);
   h2st->SetMarkerStyle(21);
   h2st->SetMarkerColor(kBlue);
   hs->Add(h2st);
   TH1D *h3st = new TH1D("h3st","test hstack",100,-4,4);
   h3st->FillRandom("gaus",10000);
   h3st->SetFillColor(kGreen);
   h3st->SetMarkerStyle(21);
   h3st->SetMarkerColor(kGreen);
   hs->Add(h3st);
   TH1D *hdata = new TH1D("hdata","test hstack",100,-4,4);
   hdata->FillRandom("gaus",10000+15000+20000);
   hdata->SetMarkerStyle(20);
   TCanvas *cst = new TCanvas("cst","stacked hists",10,10,700,700);
   hdata->Draw("e");
   hs->Draw("same");
   //   hs->Draw();
   //hdata->Draw("e,same");
   TLegend* legend = new TLegend(0.1,0.7,0.3,0.9);
   legend->AddEntry(hdata,"Data","P");
   legend->AddEntry(h1st,"first contrib.","f");
   legend->AddEntry(h2st,"2nd contrib.","f");
   legend->AddEntry(h3st,"3rd contrib.","f");
   legend->Draw();
   // dummy fit hdata->Fit("gaus");
   FitGausPlusLinear(hdata, -3., 3., 0.);
   hdata->Draw("e,same");
   return cst;
}
