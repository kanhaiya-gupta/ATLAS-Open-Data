Double_t fitfunc(Double_t* t, Double_t* par){
   Double_t x = t[0];
   Double_t sig = par[0] * TMath::Gaus(x,par[1],
                                       TMath::Abs(par[2]),kTRUE);
   Double_t bgd =par[3]* (1.+x*par[4]);
   return sig+bgd;
}

void FitGausPlusLinear(TH1D* h1, double xmin,double xmax, double xpeak=-1.){
   //gStyle->SetOptFit(111111);
   if (xpeak < 0.) xpeak=0.5*(xmin+xmax);
   TF1* fit = new TF1(TString("fit")+h1->GetName(),fitfunc,xmin,xmax, 5);
   //fit->SetLineWidth(2);
   fit->SetLineColor(h1->GetMarkerColor());
   fit->SetParNames("scale","m","#sigma","bg0","bg1");
   fit->SetParameter(0, h1->GetBinContent(h1->FindBin(xpeak)));
   fit->SetParameter(1,xpeak);
   fit->SetParameter(2,0.05*(xmax-xmin));
   fit->SetParameter(3,0.5*(h1->GetBinContent(2)+h1->GetBinContent(h1->GetNbinsX()-2)));
   h1->Fit(fit,"ER");
   TF1* fitbg=new TF1(*fit);
   fitbg->SetParameter(0,0.); 
   fitbg->SetLineColor(46);
   fitbg->Draw("same");
   TF1* fitsig=new TF1(*fit);   
   for (int i=0;i<2;i++) fitsig->SetParameter(3+i,0.); 

   h1->GetXaxis()->SetTitle("M [GeV]");
   h1->Draw("esame");
}
