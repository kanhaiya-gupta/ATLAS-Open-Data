Double_t GetBayesPrior(Double_t mu) {
  // prior probability, used below
  Double_t p=1.0;
  return p;
}

Double_t GetBayesCL(Double_t mu0,Int_t nObs,Double_t *dp=0,TRandom *rnd=0,Int_t nTry=200000) {
  /* calculates the Bayesian probability (CL) for my<=my0
     given nObs. If dp!=0 also return the error */
  Double_t sum[2],sum2[2];
  TRandom *ownedRnd=0;
  if(!rnd) ownedRnd=new TRandom3();
  if(!rnd) rnd=ownedRnd;
  sum[0]=0.; sum[1]=0.; sum2[0]=0.; sum2[1]=0.;
  /* the integral is split into two parts
     0 .. nObs               exponential distribution rising with mu
     nObs .. infinity        exponential distribution falling with mu */
  Double_t lgamNobs=TMath::LnGamma(nObs+1.);
  for(Int_t iTry=0;iTry<nTry;iTry++) {
      /* mu is drawn from an exponential distribution
	 with decay constant 1+nObs
         P(mu) = exp( -/+ (mu-nObs)/(nObs+1))
       the weight w compensates for the difference between
       the exponential distribution and the likelihood
         L(mu) = exp(-mu +ln(mu)*nObs -LnGamma(nObs+1)) 
      */
    Double_t w,mu;
    if(rnd->Uniform()>0.5) {
      mu=nObs+rnd->Exp(1.+nObs);
      w=TMath::Exp(TMath::Log(mu)*nObs-mu-lgamNobs+ (mu-nObs)/(1.+nObs));
    } else {
      mu=nObs-rnd->Exp(1.+nObs);
      if(mu>0.0) {
	w=TMath::Exp(TMath::Log(mu)*nObs-mu-lgamNobs +(nObs-mu)/(1.+nObs));
      } else {
	w=0.0;
      }
    }
    // prior probability
    w *= GetBayesPrior(mu);
    Int_t i= (mu>=mu0) ? 1 : 0;
    sum[i] += w;
    sum2[i] += w*w;
  }
  Double_t norm=sum[0]+sum[1];
  Double_t p=sum[0]/norm;
  if(dp) *dp=TMath::Sqrt(sum2[1]*sum[0]*sum[0]+sum2[0]*sum[1]*sum[1])/norm/norm;
  if(ownedRnd) delete ownedRnd;
  return p;
}

