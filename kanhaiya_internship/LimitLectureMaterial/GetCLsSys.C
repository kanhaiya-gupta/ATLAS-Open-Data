Double_t GetCLsSys(Double_t signal,Double_t bgr,
		   Double_t dBgr,Double_t dLumi,Int_t nobs,
		   Double_t *dCLs=0,TRandom *rnd=0) {
  /* calculate CLs for the given signal, background, errors
     and observed number of events, using Monte Carlo methods */
  TRandom *ownedRnd=0;
  if(!rnd) ownedRnd=new TRandom3();
  if(!rnd) rnd=ownedRnd;
  Int_t nTry=30000;
  Double_t nexp_sb=0.0;
  Double_t nexp_b=0.0;
  for(Int_t iTry=0;iTry<nTry;iTry++) {
    /* get Luminosity from truncated Gaussian */
    Double_t l=1.0;
    if(dLumi>0.0) {
      do {
	l=rnd->Gaus(1.0,dLumi);
      } while(l<=0.0);
    }
    /* get Background from truncated Gaussian */
    Double_t b=bgr;
    if(dBgr>0.0) {
      do {
	b=rnd->Gaus(bgr,dBgr);
      } while(b<=0.0);
    }
    Int_t n_b=rnd->Poisson(l*b);
    Int_t n_sb=rnd->Poisson(l*(signal+b));
    if(n_b<=nobs) nexp_b += 1.0;
    if(n_sb<=nobs) nexp_sb += 1.0;
  }
  if(ownedRnd) delete ownedRnd;
  Double_t cl_s=nexp_sb/nexp_b;
  Double_t dcl_s= TMath::Sqrt(nexp_sb+nexp_b*cl_s*cl_s)/nexp_b;
  if(dCLs) *dCLs = dcl_s;
  else {
    cout<<"CLSsys="<<cl_s<<" +/- "<<dcl_s<<" for B="<<bgr<<" +/- "<<dBgr
	<<", L=1 +/- "<<dLumi<<", signal="<<signal<<"\n";
  }
  return cl_s;
}
