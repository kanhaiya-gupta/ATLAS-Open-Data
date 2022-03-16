//#include "TRandom3.h"
//#include "TMath.h"
//using TMath::Max;
//using TMath::Min;

void fractionfitweighted()
{	
	// Example of TFractionFitter class usage
	// 1 Dimension only, x is an angle fron 0 to pi
	//
	// uses weighted histograms
	// EVT Sept. 2020
	
	// pointers to the data
	TH1D *data= new TH1D("data", "data", 10, 0., 1.);
	data->SetMarkerStyle(20);
	data->SetMarkerSize(.7);
	data->SetMinimum(0);
	TH1D *mc_wei1=new TH1D("mc1", "data", 10, 0., 1.); 
	TH1D *mc_wei2=new TH1D("mc2", "data", 10, 0., 1.);                              //data histogram, two weighted mc distributions
	TH1D *mc_unwei1=new TH1D("mc1_unwei", "data", 10, 0., 1.);   
	TH1D *mc_unwei2=new TH1D("mc2_unwei", "data", 10, 0., 1.);                               
	TH1D *weight1=new TH1D("weight1", "data", 10, 0., 1.); 
	TH1D *weight2=new TH1D("weight2", "data", 10, 0., 1.);  // weight histograms = mc_wei/mc_unwei
	mc_wei1->SetLineColor(2);
	mc_wei2->SetLineColor(3);

	int ndata=247; int nmc=2000;
	TRandom3 r;
	for (int i=0;i<ndata+nmc;i++){
		bool isdata=(i<ndata);
		bool ismc1=(r.Rndm()<0.6);
		double x=0,wei=1.;
		if (ismc1) { 
			x=r.Rndm()*r.Rndm();
			wei = 0.001+0.999; //*r.Rndm();
		}
		else {   
			x= TMath::Max(0.,TMath::Min(1.,0.5+r.Gaus(0,0.3)));
			wei= 0.5; //*(1.+r.Rndm());
		}
		if (isdata and wei>r.Rndm()) data->Fill(x);
		if (!isdata){
			if (ismc1){
				mc_wei1->Fill(x,wei);
				mc_unwei1->Fill(x);
				weight1->Fill(x,wei);
			}
			else{
				mc_wei2->Fill(x,wei);
				mc_unwei2->Fill(x);
				weight2->Fill(x,wei);
			}
		}
	}
	
	weight1->Divide(mc_unwei1);
	weight2->Divide(mc_unwei2);
	// FractionFitter
	TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
	mc->Add(mc_wei1);
	mc->Add(mc_wei2);
	TFractionFitter* fit = new TFractionFitter(data, mc); // initialise
	fit->Constrain(0,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	fit->Constrain(1,0.0,1.0);               // constrain fraction 1 to be between 0 and 1
	fit->SetWeight(0,weight1);
	fit->SetWeight(1,weight2);
	//fit->SetRangeX(1,15);                    // use only the first 15 bins in the fit
	Int_t status = fit->Fit();               // perform the fit
	cout << "fit status: " << status << endl;
	

	// Display
	gStyle->SetOptStat(0);
	TCanvas* c1= new TCanvas("c1", "FractionFitter example", 700, 700);
	c1->Divide(2,2);
	c1->cd(1);
	mc_wei1->Draw();
	mc_wei2->Draw("same");
	c1->cd(2);
	weight1->Draw();
	weight2->Draw("same");
	c1->cd(3);
	data->Draw("e");
        c1->cd(4);
	if (status == 0) {
	  TH1F* result = (TH1F*) fit->GetPlot();
	  data->Draw("Ep");
	  result->Draw("same");}

	
}
