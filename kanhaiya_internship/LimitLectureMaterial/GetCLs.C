Double_t GetCLs(Double_t signal,Double_t bgr,Int_t nobs) {
  /* calculate CLs for the given
     signal, background , nobs. Does not work for high nobs. */
  Double_t cl_s;
  Double_t cl_b=0.0;
  Double_t cl_sb=0.0;
  Double_t lnGamma=0.0;
  Double_t logB=TMath::Log(bgr);
  Double_t logSB=TMath::Log(signal+bgr);
  for(Int_t i=0;i<=nobs;i++) {
    Double_t p_b=TMath::Exp(i*logB-bgr-lnGamma);
    Double_t p_sb=TMath::Exp(i*logSB-signal-bgr-lnGamma);
    cl_b += p_b;
    cl_sb += p_sb;
    lnGamma +=TMath::Log(i+1);
  }
  cl_s=cl_sb/cl_b;
  return cl_s;
}
