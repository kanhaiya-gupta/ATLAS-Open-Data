#include <cstdlib> 
#include <cstdio> 
#include <string> 
#include <iostream> 
#include <fstream>
#include "TSystem.h"
#include "myroot.h"
TH1D* StackPlot(TString dir, TString name){

   TFile *fileData = TFile::Open("data.root");
   const int nmcfiles=3;
   TString filenames[nmcfiles]={"ZZ.root","WZ.root","BSM1.root"};
   int colors[nmcfiles]={kGreen, kOrange, kRed};
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
         leg->AddEntry(hist,filenames[i],"F");
      }
   }
   double max=hdata->GetMaximum();
   if (mc->GetMaximum()>hdata->GetMaximum()) max=mc->GetMaximum();
   hdata->SetMaximum(1.3*max);
   hdata->Draw("e");   
   mc->Draw("samehisto"); 
   hdata->Draw("esame"); 
   leg->Draw();
   return hdata;
}
