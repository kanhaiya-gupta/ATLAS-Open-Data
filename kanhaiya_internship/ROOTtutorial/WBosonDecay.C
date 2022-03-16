// usage
// root -l
// .L WBosonDecay.C
//  WBosonDecay()->Draw()
TH1D* WBosonDecay(){
  TRandom3 r;
  TH1D* h = new TH1D("ptlep","ptlep",100,0.,100.);
  TLorentzVector pWr, plepr,pnur; // vectors in W rest frame;
    TLorentzVector pW, plep,pnu;
  double mW=80.3;
  double twopi = 2.*acos(-1.);

  for (int itry=0;itry<10000;itry++){
    double pz = r.Gaus(0.,50.);
    mW = r.BreitWigner(80.3,2.); // finite width of W boson
    pW.SetPxPyPzE(0,0,pz,sqrt(pz*pz+mW*mW));
    TVector3 boost = pW.BoostVector();
    double thet=acos(2.*r.Rndm()-1.);
    double phi =twopi * r.Rndm();
    TVector3 mom; mom.SetMagThetaPhi(mW/2., thet,phi);
    plep.SetVectM(mom,0);pnu.SetVectM(-mom,0);
    if (plep.M()>0.001) cout << "Error m="<<plep.M()<<endl;
    plep.Boost(boost);
    pnu.Boost(boost);
    h->Fill(plep.Pt());
  }
  return h;
}
