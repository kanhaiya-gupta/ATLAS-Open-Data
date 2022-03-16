#include "GetPoissonLimit.C"

Double_t GetExpectedLimit(Double_t cl,Double_t mu) {
  /* calculate the expected limit, given the confidence level
     and the Poisson parameter mu. Note: does not work for very large mu */
  Int_t n0=(int)mu;
  Int_t n1=n0+1;
  Double_t p0=1.0;
  if(mu>0.0) p0=TMath::Exp(n0*TMath::Log(mu)-mu-TMath::LnGamma(n0+1.));
  Double_t p1=p0*mu/n1;
  Double_t sum=0.0;
  while((p0>0.)&&(n0>=0)) {
    sum += p0*GetPoissonLimit(cl,n0);
    p0 *= n0/mu;
    n0--;
  }
  while(p1>0.0) {
    sum += p1*GetPoissonLimit(cl,n1);
    n1++;
    p1 *= mu/n1;
  }
  return sum;
}
