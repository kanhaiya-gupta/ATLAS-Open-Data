Double_t GetPoissonLimit(Double_t cl,Int_t n0) {
  /* calculates the limit on the number of events
     given n0 */
  Double_t mu=0.5*TMath::ChisquareQuantile(cl,2*(n0+1));
  return mu;
}
