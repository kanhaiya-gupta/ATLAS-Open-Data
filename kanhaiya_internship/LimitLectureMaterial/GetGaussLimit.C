Double_t GetGaussLimit(Double_t cl,Int_t n0) {
  /* calculates the limit on the number of events 
     given n0 in the Gaussian approximation. */
  Double_t mu=TMath::ErfcInverse(2.*(1.-cl))*TMath::Sqrt(2.*n0)+n0;
  return mu;
}
