#include <cstdlib> 
#include <cstdio> 
#include <string> 
#include <iostream> 
#include <fstream>
#include "TSystem.h"
#include "myroot.h"
TH1D* StackPlot(TString dir, TString name){

   gStyle->SetOptStat(0);
   TFile *fileData = TFile::Open("data.root");
   const int nmcfiles=4;
   TString filenames[nmcfiles]={"WZ1.root","Zjet1.root","ZZ1.root","ttbar1.root"};
   TString file_names[nmcfiles]={"WZ","Z+jet","ZZ","t#bar{t}"};
   int colors[nmcfiles]={kOrange,kGray,kGreen,kRed};
   TFile* files[nmcfiles];
   THStack *mc = new THStack("mc","MC background");
   TLegend* leg = new TLegend(0.15,0.65,0.35,0.89);
   TString fullname;
   if (dir =="") fullname = name;
   else fullname = dir + "/" + name;
   TH1D* hdata = (TH1D*) fileData->Get(fullname);
   if (!hdata) return 0;
   hdata->SetMarkerStyle(20);
   hdata->SetMarkerSize(1.);
   hdata->SetMarkerColor(1);
   //  hdata->GetXaxis()->SetTitle("P_{T}^{lead} (GeV/c)");
   leg->AddEntry(hdata, "Data", "PE");
   for (int i=0;i<nmcfiles;i++){
      TFile* fil = (TFile*) gROOT->FindObject(filenames[i]);
      if (!fil){
         fil = TFile::Open(filenames[i]);
      }
      TH1D* hist = (TH1D*) fil->Get(fullname);
      if (hist){
         hist->SetFillColor(colors[i]);
         mc->Add(hist);
         leg->AddEntry(hist,file_names[i],"F");
      }
   }
   double max=hdata->GetMaximum();
   if (mc->GetMaximum()>hdata->GetMaximum()) max=mc->GetMaximum();
   hdata->SetMaximum(1.3*max);
   hdata->Draw("e");   
   mc->Draw("samehisto"); 
   hdata->Draw("esame");
   leg->SetLineColor(0); 
   leg->Draw();
   return hdata;
}
