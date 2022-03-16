Double_t GetBayesBgrPrior(Double_t sig) {
  // prior probability, used below
  Double_t p=(sig>=0.) ? 1. : 0.;
  return p;
}

Double_t GetBayesBgrCL(Double_t sig0,Double_t bgr,Int_t n0,Double_t *dp=0,
		       TRandom *rnd=0) {
  /* calculates the Bayesian probability (1-CL) for my>=sig0+bgr
     given n0. If dp!=0 also return the error */
  Int_t nTry=100000;
  Double_t sum[2],sum2[2];
  TRandom *ownedRnd=0;
  if(!rnd) ownedRnd=new TRandom3();
  if(!rnd) rnd=ownedRnd;
  sum[0]=0.; sum[1]=0.; sum2[0]=0.; sum2[1]=0.;
  /* the integral is split into two parts
     0 .. n0               exponential distribution rising with mu
     n0 .. infinity        exponential distribution falling with mu */
  Double_t lgamN0=TMath::LnGamma(n0+1.);
  for(Int_t iTry=0;iTry<nTry;iTry++) {
      /* mu is drawn from an exponential distribution
	 with decay constant 1+n0
         P(mu) = exp( -/+ (mu-n0)/(n0+1))
       the weight w compensates for the differnce between
       the exponential distribution and the likelihood
         L(mu) = exp(-mu +ln(mu)*n0 -LnGamma(n0+1)) 
      */
    Double_t w,mu;
    if(rnd->Uniform()>0.5) {
      mu=n0+rnd->Exp(1.+n0);
      w=TMath::Exp(TMath::Log(mu)*n0-mu-lgamN0+ (mu-n0)/(1.+n0));
    } else {
      mu=n0-rnd->Exp(1.+n0);
      w=TMath::Exp(TMath::Log(mu)*n0-mu-lgamN0 +(n0-mu)/(1.+n0));
    }
    // reject mu smaller zero
    if(mu<0.0) continue;
    Double_t sig=mu-bgr;
    // prior probability
    w *= GetBayesBgrPrior(sig);
    Int_t i= (sig>=sig0) ? 1 : 0;
    sum[i] += w;
    sum2[i] += w*w;
  }
  Double_t norm=sum[0]+sum[1];
  Double_t p=sum[0]/norm;
  if(dp) *dp=TMath::Sqrt(sum2[1]*sum[0]*sum[0]+sum2[0]*sum[1]*sum[1])/norm/norm;
  if(ownedRnd) delete ownedRnd;
  return p;
}
