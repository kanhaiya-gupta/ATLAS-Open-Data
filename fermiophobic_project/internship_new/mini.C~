#define mini_cxx
#include "mini.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
using std::vector;
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"

void mini::Loop()
{
//   In a ROOT session, you can do:
//      root> .L mini.C
//      root> mini t
//
// full processing of all ATLAS 2lepton data in a root session
// root>  TChain ch("mini")
//   ch.Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/Data/data_*.2lep.root")
//   .L mini.C
//   mini t(ch)
//   t.Loop()
//   
//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
   if (fChain == 0) return;

   //std::ofstream myfile;
   //myfile.open("2lep.csv");
   //std::cout << " created 2lep.csv" <<std::endl;
   Long64_t nentries = fChain->GetEntriesFast();
   // nentries = 1000000;  // number of events
   Long64_t mu3count = 0;
   Long64_t zmm_count = 0;
   Long64_t nbytes = 0, nb = 0;
   float mw = 80.4;
  
  
   // booking histograms
   TFile outf("output.root","RECREATE");  // create a root file
   std::cout << "creating output.root" <<std::endl;
   
   // Defining Histograms
   TH1D* h_mupt= new TH1D("mupt","Muon Transverse Momentum; P_{T} (GeV/c); Events per bin",200,0.,100.);
   TH1D* h_mdimu= new TH1D("mass_dimu","Di-muon invariant mass; M_{#mu^{+}#mu^{-}}(GeV/c^{2});Events per bin",400,0.,200.);
   TH1D* h_mtrimu= new TH1D("mass_trimu","Tri-muon invariant mass; M (GeV/c^{2});Events per bin",200,0.,2000.);
   TH1D* h_deltaR= new TH1D("deltaR","deltaR;#Delta R; Events per bin",400,0.,10.0);
   TH1D* h_charge= new TH1D("sum_charge","sum_charge; #Sigma_{q_{i}};Events per bin",40,-10,10);
   TH1D* h_ptcone1 = new TH1D("ptconemu1","ptconemu1", 400,0.,0.8);
   TH1D* h_ptcone2 = new TH1D("ptconemu2","ptconemu2", 400,0.,0.8);
   TH1D* h_ptcone3 = new TH1D("ptconemu3","ptconemu3", 400,0.,0.8);
   TH1D* h_met = new TH1D("met","met; Missing Transverse Energy; #slash{E_{T}} (GeV); Events per bin",100,0.,600.);
   TH1D* h_Mt = new TH1D("M_T(W)","Transverse mass ;M_{T}[W] (GeV/c^{2});Events per bin",50,30.,1000.);
   TH1D* h_Mt_wz = new TH1D("M_T(WZ)","Transverse mass of WZ; M_{T}[WZ](GeV/c^{2}); Events per bin",80,40.,1000.);
   TH1D* h_M_wz = new TH1D("M_(WZ)","Invariant mass of WZ; M_{inv}[WZ](GeV/c^{2});Events per bin",80,40.,1000.);
   TH1D* h_nan = new TH1D("N_nan","Histogram plot of  the argument of cosh^-1(..) BSM3 WZ_500p;;Events per bin", 80,-10., 30.);
   TH2D* h2_plot1 = new TH2D("2d_plot1","Et_mis vs M_T(WZ);M_T(WZ);Et_mis",50,0.,800.,50, 10., 200.);
   TH2D* h2_plot2 = new TH2D("2d_plot2","",50,0.,800.,50,10.,200.);

   TH2D* h_tri = new TH2D("2d_tri","2D Et_mis vs Tri-muons inv mass;M_{3#muon}(inv);Et_mis",100,0.,800.,100,0.,800.);

   TH1D* h_phi = new TH1D("d_phi","Cos(#Delta #Phi) plot of #mu #nu;;Events per bin",80,-5.,5.);
   TH1D* h_eta = new TH1D("d_eta","#Delta #eta plot;;Events per bin",80,-10.,10.);
   TH1D* h_arg = new TH1D("d_arg","m_w^2/(2*Pt*Etmis);;Events per bin",80,-2.,20.);


    TH1D* h_Mt_weit = new TH1D("M_T_weit(W)","Transverse mass ;M_{T}[W] (GeV/c^{2});Events per bin",25,35.,120.);
    TH1D* h_Mt_wz_weit = new TH1D("M_T_weit(WZ)","Transverse mass of WZ; M_{T}[WZ](GeV/c^{2}); Events per bin",80,40.,1000.);

    TH1D* h_Mt_weit_mc = new TH1D("M_T_weit_mc(W)","Transverse mass ;M_{T}[W] (GeV/c^{2});Events per bin",25,35.,120.);
    TH1D* h_Mt_wz_weit_mc = new TH1D("M_T_weit_mc(WZ)","Transverse mass of WZ; M_{T}[WZ](GeV/c^{2}); Events per bin",80,40.,1000.);

    TH1D* h_deltaRnu = new TH1D("Delta_Rnu","#Delta#Phi between #mu & #nu;#Delta#Phi; Events per bin",80,-10,10.);

    TH1D* h_mtrimu_bkg = new TH1D("mass_trimu_bkg","Tri-muon invariant mass; M (GeV/c^{2});Events per bin",200,0.,2000.);
    TH1D* h_mtrimu_sig = new TH1D("mass_trimu_sig","Tri-muon invariant mass; M (GeV/c^{2});Events per bin",200,0.,2000.);
   
   // di-mons histograms
   /*  
   TH1D* d1_mupt= new TH1D("d1-mupt","d1-mupt",200,0.,100.);
   TH1D* d2_mupt= new TH1D("d2-mupt","d2-mupt",200,0.,100.);
   TH1D* d_mdimu= new TH1D("d-mass_dimu","d-mass_dimu",400,0.,200.);
   TH1D* d1_Mt = new TH1D("d1-M_T","d1-M_T",400,0.,200.);
   TH1D* d2_Mt = new TH1D("d2-M_T","d2-M_T",400,0.,200.);
   TH1D* dz_Mt = new TH1D("dz-M_T","dz-M_T",400,0.,200.);
   */ 
   
   for  (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // calculate M?T for the firSt muon in di-lepton events
      /* 
      if (lep_n != 2) continue;
      if ((*lep_type)[0] != 11) continue; // select 3 mu
      if ((*lep_type)[1] != 13) continue; // select 3 mu
      if ((*lep_pt)[0]>13.e6 or (*lep_pt)[1]>13.e6) continue; // discard mismeasured events
      if ((*lep_isTightID)[0] == false) continue; 
      if ((*lep_isTightID)[1] == false) continue;
       if ((*lep_pt)[0]<20000.) continue;
      if ((*lep_pt)[1]<20000.) continue;

      TLorentzVector dmu1,dmu2,MeT1; // 4-momentum in GeV
	    
      dmu1.SetPtEtaPhiM((*lep_pt)[0]/1000.,(*lep_eta)[0],(*lep_phi)[0],0.1057);
      dmu2.SetPtEtaPhiM((*lep_pt)[1]/1000.,(*lep_eta)[1],(*lep_phi)[1],0.0005);
      MeT1.SetPtEtaPhiE((met_et)/1000.,0.,met_phi,(met_et)/1000.);

      float dmu12 = (dmu1 + dmu2).M();
       
       d1_mupt->Fill(dmu1.Pt());
       d2_mupt->Fill(dmu2.Pt());
       d_mdimu->Fill(dmu12);

       float mtw1 = 0;
       float mtw2 = 0;
       float mtz = 0;

       if (TMath::Abs(dmu12 - 92.1) > 15)
      {  mtw1 = sqrt(2*dmu1.Pt()*MeT1.E()*(1-cos(dmu1.DeltaPhi(MeT1))));
       mtw2 = sqrt(2*dmu2.Pt()*MeT1.E()*(1-cos(dmu2.DeltaPhi(MeT1))));
       mtz = sqrt(2*dmu1.Pt()*dmu2.Pt()*(1-cos(dmu1.DeltaPhi(dmu2))));
       
       
       d1_Mt->Fill(mtw1);
       d2_Mt->Fill(mtw2);
       dz_Mt->Fill(mtz);
      }

       // DI MUON END

       */
     
      
      				 
      // select e mu final state
      if (lep_n != 3) continue;
      if ((*lep_type)[0] != 13 or (*lep_type)[0] != 11) continue; // select 3 mu
      if ((*lep_type)[1] != 13) continue; // select 3 mu
      if ((*lep_type)[2] != 13) continue; // select 3 mu
      if ((*lep_pt)[0]>13.e6 or (*lep_pt)[1]>13.e6 or (*lep_pt)[2]>13.e6) continue; // discard mismeasured events

      // tight lepton selection
      if ((*lep_isTightID)[0] == false) continue; 
      if ((*lep_isTightID)[1] == false) continue;
      if ((*lep_isTightID)[2] == false) continue;

      // hard leptons
      if ((*lep_pt)[0]<20000.) continue;
      if ((*lep_pt)[1]<20000.) continue;
      if ((*lep_pt)[2]<20000.) continue;
      
      TLorentzVector pmu1,pmu2,pmu3,MeT,pnu1, pnu2; // 4-momentum in GeV
      pmu1.SetPtEtaPhiM((*lep_pt)[0]/1000.,(*lep_eta)[0],(*lep_phi)[0],0.1057);
      pmu2.SetPtEtaPhiM((*lep_pt)[1]/1000.,(*lep_eta)[1],(*lep_phi)[1],0.1057);
      pmu3.SetPtEtaPhiM((*lep_pt)[2]/1000.,(*lep_eta)[2],(*lep_phi)[2],0.1057);
      MeT.SetPtEtaPhiM((met_et)/1000.,0.,met_phi,0.);

       for(unsigned int i=0; i<lep_n; i++)
		    {
                      TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);
		    }
      
      float m2mu = 0;
      float m3mu = (pmu1+pmu2+pmu3).M();
      float Delta_m1m2 = pmu1.DeltaR(pmu2);
      float Delta_m1m3 = pmu1.DeltaR(pmu3);
      float Delta_m2m3 = pmu2.DeltaR(pmu3);
      float max_R = TMath::Max(Delta_m1m2, Delta_m1m3);
      float max_delR = TMath::Max(max_R, Delta_m2m3);
      int lep_char = (*lep_charge)[0] + (*lep_charge)[1] + (*lep_charge)[2];
      float pt_cone1 = (*lep_ptcone30)[0]/(*lep_pt)[0];
      float pt_cone2 = (*lep_ptcone30)[1]/(*lep_pt)[1];
      float pt_cone3 = (*lep_ptcone30)[2]/(*lep_pt)[2];

      //Scale factors
      // luminosity = N/XSection = nentries*mc/XSection;
      // To make the luminosity same, sacle = (lumi of data)/(lumi of MC);
      
      float lumi_mc = 10.*XSection*1000./(nentries*mcWeight);
      float scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP*XSection/1000;
      
      //MC weight
      Float_t m_mcWeight = mcWeight;

      //Total weight
      //   float weight = 2.5*scaleFactor*m_mcWeight;    // XSection/1000.;

     float weight = 1. ; // XSection/1000.*m_mcWeight; // ZZ   // for data  0.055*m_mcWeight; WZ
      // Transverse mss of W boson Calculation

      float m12 = (pmu1 + pmu2).M();
      float m13 = (pmu1 + pmu3).M();
      float m23 = (pmu2 + pmu3).M();

      float delta_m12 =-1.;
      float delta_m13 =-1.;
      float delta_m23 =-1.;
     

      // leptons coming from the z have different charge

      if ( (*lep_charge)[0]* (*lep_charge)[1] < 0)
	{ delta_m12 = TMath::Abs(m12-91.2);}

       if ( (*lep_charge)[0]* (*lep_charge)[2] < 0)
	{ delta_m13 = TMath::Abs(m13-91.2);}

        if ( (*lep_charge)[1]* (*lep_charge)[2] < 0)
	{ delta_m23 = TMath::Abs(m23-91.2);}

      // define candidates 
	 int wcand = 0;
	 float tmp = 0.;
	

	  // in the case we have eee or mumumu
	      ////		      
	      if( ( delta_m12 >0 && delta_m23 >0)  && delta_m13 < 0 && (delta_m12 < delta_m23) ) {
		tmp = delta_m12; 
		wcand = 3;
		zmm_count++;
	      } 
	      
	      if( ( delta_m12 >0 && delta_m23 >0)   && delta_m13 < 0 && (delta_m12 > delta_m23) ) {
		tmp = delta_m23; 
		wcand = 1;
		zmm_count++;
	      }
	      
	      if( ( delta_m12 >0 && delta_m13 >0 )  && delta_m23 < 0 && (delta_m12 < delta_m13) ) {
		tmp = delta_m12;
		wcand = 3;
		zmm_count++;
	      }
	      
	      if( ( delta_m12 >0 && delta_m13 >0 ) && delta_m23 < 0 && (delta_m12 > delta_m13) ) {
		tmp = delta_m13;
		wcand = 2;
		zmm_count++;
	      }
	      
	      if( ( delta_m13 >0 && delta_m23 >0)  && delta_m12 < 0 && (delta_m13 < delta_m23) ) {
		tmp = delta_m13;
		wcand = 2;
		zmm_count++;
	      }
	      
	      if( ( delta_m13 >0 && delta_m23 >0)  && delta_m12 < 0 && (delta_m13 > delta_m23) ) {
		tmp = delta_m23;
		wcand = 1;
		zmm_count++;
	      }
	 

      float mtw = 0.;
      float eta_1 = 0.;
      float eta_2 = 0.;
      float mwz = 0.;
      float arg_cosh = 0.;
      float d_phi = 0;
      float d_eta = 0.;
      float arg_1 = 0.;

       if  (wcand == 1) 
       	{  mtw = sqrt(2*pmu1.Pt()*MeT.Pt()*(1-cos(pmu1.DeltaPhi(MeT))));
	   m2mu =  (pmu2+pmu3).M();
	  eta_1 = pmu1.Eta() - TMath::ACosH((mw*mw)/(2*pmu1.Pt()*MeT.E()) + cos(pmu1.DeltaPhi(MeT)));
	  eta_2 = pmu1.Eta() + TMath::ACosH((mw*mw)/(2*pmu1.Pt()*MeT.E()) + cos(pmu1.DeltaPhi(MeT)));
	  pnu1.SetPtEtaPhiM((met_et)/1000.,eta_1,met_phi,0.);
	  pnu2.SetPtEtaPhiM((met_et)/1000.,eta_2,met_phi,0.);
	  mwz = (pmu1 + pmu2 +  pmu3 + pnu1).M();
          arg_cosh =  cos(pmu1.DeltaPhi(MeT)); // (mw*mw)/(2*pmu1.Pt()*MeT.E()) + cos(pmu1.DeltaPhi(MeT));

	  d_phi = cos(pmu1.DeltaPhi(MeT));
	  d_eta = pmu1.Eta() - MeT.Eta();
	  arg_1 = (mw*mw)/(2*pmu1.Pt()*MeT.E());

	   h_deltaRnu->Fill(pmu1.DeltaPhi(MeT));
	  
	  //  mwz = (pmu1 + pmu2 +  pmu3 + pnu2).M();
	   
	}

        if  (wcand == 2)
       	{  mtw = sqrt(2*pmu2.Pt()*MeT.Pt()*(1-cos(pmu2.DeltaPhi(MeT))));
	   m2mu =  (pmu1+pmu3).M();
	   eta_1 = pmu2.Eta() - TMath::ACosH((mw*mw)/(2*pmu2.Pt()*MeT.E()) + cos(pmu2.DeltaPhi(MeT)));
           eta_2 = pmu2.Eta() + TMath::ACosH((mw*mw)/(2*pmu2.Pt()*MeT.E()) + cos(pmu2.DeltaPhi(MeT)));
	   pnu1.SetPtEtaPhiM((met_et)/1000.,eta_1,met_phi,0.);
	   pnu2.SetPtEtaPhiM((met_et)/1000.,eta_2,met_phi,0.);
	   mwz = (pmu1 + pmu2 +  pmu3 + pnu1).M();
	   arg_cosh =  cos(pmu2.DeltaPhi(MeT)); // (mw*mw)/(2*pmu1.Pt()*MeT.E()) + cos(pmu2.DeltaPhi(MeT));

	   d_phi = cos(pmu2.DeltaPhi(MeT));
	   d_eta = pmu2.Eta() - MeT.Eta();
	   arg_1 = (mw*mw)/(2*pmu2.Pt()*MeT.E());
	   //   mwz = (pmu1 + pmu2 +  pmu3 + pnu2).M();
	   h_deltaRnu->Fill(pmu2.DeltaPhi(MeT));
	}
       
      if  (wcand == 3)
	 {  mtw = sqrt(2*pmu3.Pt()*MeT.Pt()*(1-cos(pmu3.DeltaPhi(MeT))));
	    m2mu =  (pmu1+pmu2).M();
            eta_1 = pmu3.Eta() - TMath::ACosH((mw*mw)/(2*pmu3.Pt()*MeT.E()) + cos(pmu3.DeltaPhi(MeT)));
	    eta_2 = pmu3.Eta() + TMath::ACosH((mw*mw)/(2*pmu3.Pt()*MeT.E()) + cos(pmu3.DeltaPhi(MeT)));
	    pnu1.SetPtEtaPhiM((met_et)/1000.,eta_1,met_phi,0.);
	    pnu2.SetPtEtaPhiM((met_et)/1000.,eta_2,met_phi,0.);
	    mwz = (pmu1 + pmu2 +  pmu3 + pnu1).M();
	    arg_cosh = cos(pmu1.DeltaPhi(MeT)); //  (mw*mw)/(2*pmu1.Pt()*MeT.E()) + cos(pmu1.DeltaPhi(MeT));

	    d_phi = cos(pmu1.DeltaPhi(MeT));
	    d_eta = pmu3.Eta() - MeT.Eta();
	    arg_1 = (mw*mw)/(2*pmu3.Pt()*MeT.E());
	    //  mwz = (pmu1 + pmu2 +  pmu3 + pnu2).M();

	    h_deltaRnu->Fill(pmu3.DeltaPhi(MeT));
	 }


      // cut: m_ll - m(Z) < 10
	      if(tmp < 10.)
		{
		  // mtw > 30 GeV, at least one lepton with pT > 25 GeV 
		  if( mtw > 30. && MeT.Pt() > 30. && (pmu1.Pt() > 25. || pmu2.Pt() > 25. || pmu3.Pt() > 25.) )
		    {
		      h_Mt->Fill(mtw, weight);
		      h_Mt_weit->Fill(mtw,weight);
		     
		       h_Mt_weit_mc->Fill(mtw,weight);
		  
		      }}
	     
	     
	      h_M_wz->Fill(mwz,weight);
	      h_mupt->Fill(pmu1.Pt(),weight);
	      h_mupt->Fill(pmu2.Pt(),weight);
	      h_mupt->Fill(pmu3.Pt(),weight);
	      h_mdimu->Fill(m2mu,weight);
      // if (Delta_m1m2 < 0.4 && Delta_m1m3 < 0.4 && Delta_m2m3 < 0.4)
   
      h_mtrimu->Fill(m3mu,weight);
      h_deltaR->Fill(max_delR,weight);
      h_charge->Fill(lep_char,weight);
      h_ptcone1->Fill(pt_cone1,weight);
      h_ptcone2->Fill(pt_cone2,weight);
      h_ptcone3->Fill(pt_cone3,weight);
      h_met->Fill(MeT.Pt(),weight);
      h_nan->Fill(arg_cosh,weight);

      h_phi->Fill(d_phi,weight);
      h_eta->Fill(d_eta,weight);
      h_arg->Fill(arg_1,weight);

      h_tri->Fill(MeT.Pt(), m3mu);

      if ((MeT.Pt()>=30 && MeT.Pt()<=80) && (m3mu>=95 && m3mu<=400))
	{ h_mtrimu_bkg->Fill(m3mu,weight);}

      if ((MeT.Pt()>80 && MeT.Pt()<=500) && (m3mu>400 && m3mu<=2000))
	{ h_mtrimu_sig->Fill(m3mu,weight);}
     
      // Calculation for the Transverse mass of WZ boson pair 

      float p_x1 = pmu1.Pt()*TMath::Cos(pmu1.Phi());
      float p_x2 = pmu2.Pt()*TMath::Cos(pmu2.Phi());
      float p_x3 = pmu3.Pt()*TMath::Cos(pmu3.Phi());

      float p_y1 = pmu1.Pt()*TMath::Sin(pmu1.Phi());
      float p_y2 = pmu2.Pt()*TMath::Sin(pmu2.Phi());
      float p_y3 = pmu3.Pt()*TMath::Sin(pmu3.Phi());

      float et_x = MeT.Pt()*TMath::Cos(MeT.Phi());
      float et_y = MeT.Pt()*TMath::Sin(MeT.Phi());

      float mtwz = sqrt(pow((pmu1.Pt() + pmu2.Pt() + pmu3.Pt() + MeT.E()),2)  - pow((p_x1 + p_x2 + p_x3 + et_x ),2) -  pow((p_y1 + p_y2 + p_y3 + et_y ),2));

      if ((met_et/1000 >30.&& met_et/1000. < 80.) && (mtwz > 120 && mtwz <350)){
        h_Mt_wz->Fill(mtwz);
	h_Mt_wz_weit->Fill(mtwz,weight);
      }
      

      //  h2_plot1->Fill(mtwz,met_et/1000.);

       
      
      if (met_et/1000 >80. && mtwz > 350)
	{ h2_plot1->Fill(mtwz,met_et/1000.,weight);
	 
	}
      if (met_et/1000. <=80. && mtwz <= 350.)
	{ h2_plot2->Fill(mtwz,met_et/1000.,weight);
         
	}
       
      
      
       //  h2_plot1->Draw("same box");
      
      // Invariant mass Calcuation
     
      

      
     
      
       
      //if (mu3count==0) myfile<<"pt1,pt2,pt3,m3mu\n";
      //if (mu3count<10) std::cout<<pmu1.Pt()<<","<<pmu2.Pt()<<","<<pmu3.Pt()<<","<<m3mu<<std::endl;
      //myfile<<pmu1.Pt()<<","<<pmu2.Pt()<<","<<pmu3.Pt()<<","<<m3mu<<std::endl;
      //mu3count++;
   }

   h_Mt_weit_mc->Divide(h_Mt);
   std::cout<<"Number of Z-->mumu + 1mu events: "<<zmm_count<<std::endl;
   
   //   float h1scale = h_Mt_s->Integral();
   // float scale1 = 1/h1scale;
   // h_Mt_s->Scale(scale1);

   // float h2scale =  h_Mt_wz_s->Integral();
   //float scale2 = 1/h2scale;
   //h_Mt_wz_s->Scale(scale2);
   
   outf.Write();
   //myfile.close();
   
}

mini::mini(TTree *tree) : fChain(0) 
 {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
     // TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("data_A.2lep.root");
     // if (!f || !f->IsOpen()) {
     //    f = new TFile("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/Data/data_A.2lep.root");
     // }
     //  f->GetObject("mini",tree);
      TChain* tchain = new TChain("mini");
      //   tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/Data/data_*.2lep.root"); 
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363491.lllv.2lep.root"); // WZ  MC_actual
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363490.llll.2lep.root");  // ZZ MC
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_3641*.Zmumu_PTV0_70_CVetoBVeto.2lep.root");  //  Z + jet MC  
          tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392217.C1N2_WZ_400p0_0p0_3L_2L7.2lep.root"); // M1 
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392220.C1N2_WZ_350p0_0p0_3L_2L7.2lep.root"); // M2 
      //    tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392223.C1N2_WZ_500p0_0p0_3L_2L7.2lep.root"); // M3 
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392226.C1N2_WZ_100p0_0p0_3L_2L7.2lep.root"); // M4
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392302.C1N2_WZ_500p0_100p0_2L2J_2L7.2lep.root"); // M5
      //   tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_301325.ZPrime1000_tt.2lep.root");   // BSM Z' --> tt~
      
	 
      tree = tchain;
   }
   Init(tree);
   Loop();
}
