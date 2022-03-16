#include "myroot.h"

void histo1()
{
   cout << "simple example on how to use histograms"<<endl;
   TRandom3 r;
   TFile* outfile = new TFile("histo1.root","RECREATE");
   TH1D* h = new TH1D("example", "example",50,0.,1.);

   for (int i=0;i<1000;i++){
      double x= r.Rndm();
      double y = x*x;
      h->Fill(y);
   }

   h->Write();

   h->Draw();
}


