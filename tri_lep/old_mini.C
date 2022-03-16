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

   if (fChain == 0) return;

  
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   // nentries = 10000;  // number of events

   int count_init = 0;
   int count_trilep = 0;
   int count_met = 0;
   int count_trig = 0;
   int count_bjet = 0;
   int count_mll = 0;
   int count_mbb = 0;
   int count_mllbb = 0;
   
  
 
   // booking histograms
   TFile outf("output.root","RECREATE");  // create a root file
   std::cout << "creating output.root" <<std::endl;
   
   // Defining Histograms
   //   TH1D* h_leppt= new TH1D("leppt","Leading P_{T} lepton; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   //   TH1D* h_jetpt= new TH1D("jetpt","Leading P_{T} b-jet; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   //   TH1D* h_subjetpt= new TH1D("subjetpt","Sub-leading P_{T} b-jet; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   TH1D* h_mtw= new TH1D("mtw","Transverse mass of W boson; M_{T} (GeV/c^{2}); Events per bin",80,20.,150.);
   TH1D* h_wmass= new TH1D("wmass","Invariant mass of W boson; M_{jj} (GeV/c^{2}); Events per bin",80,0.,200.);
   TH1D* h_mll= new TH1D("zmass","Invariant mass of Z boson; M_{ll} (GeV/c^{2}); Events per bin",100,20.,150.);
   TH1D* h_pnuz= new TH1D("pnuz","P_{z} of neutrino; P_{T} (GeV/c); Events per bin", 200,-1000., 5000.);
  //  TH2D* h_lepjet2d = new TH2D("lepjet2d","2D Histogram plot ;M_{ll}/92 (GeV/c^{2});M_{bb}/125 (GeV/c^{2})",100,0.6,1.4,100, 0.2, 2.);

  
   
   for  (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      count_init = count_init + 1;

   
     //Scale factors
      Float_t scaleFactor = scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP;

      //MC weight
      float  m_mcWeight = mcWeight;

      //Total weight
      //   float weight = scaleFactor*m_mcWeight*XSection/1000.;

      float weight = 1.;    // data


  
//Preselection cut for electron/muon trigger
	if(trigE || trigM){

      count_trig = count_trig + 1; 
     
	  // Preselection of good leptons
	  int goodlep_index[3];
	  int goodlep_n = 0;
	  int lep_index =0;
	  
	  for(unsigned int i=0; i<lep_n; i++)
	    {
              TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);

	      // Lepton is Tight
	      if( lep_isTightID->at(i) )
		{
		  //  Lepton is isolated and with at least 20 GeV
		  if( lep_pt->at(i) > 20000. && ( (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.15) && ( (lep_etcone20->at(i) / lep_pt->at(i)) < 0.15) )
		    {
		      if ( lep_type->at(i)==11 && TMath::Abs(lep_eta->at(i)) < 2.47 && ( TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52 ) ) {
                      if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
                        goodlep_n = goodlep_n + 1;
                        goodlep_index[lep_index] = i;
                        lep_index++;
                      }
                 }
      // muon selection

           // muon selection
                      if ( lep_type->at(i) == 13 && TMath::Abs(lep_eta->at(i)) < 2.5 ) {
                        if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 3 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
                          goodlep_n = goodlep_n + 1;
                          goodlep_index[lep_index] = i;
                          lep_index++;
                        }
                      }
                    }
                }
            }

	  
		  
	  //Exactly three good leptons
	  if(goodlep_n==3 )
	    {
	      count_trilep =  count_trilep + 1;
	      
	      int goodlep1_index = goodlep_index[0];
	      int goodlep2_index = goodlep_index[1];
	      int goodlep3_index = goodlep_index[2];
	      
	      
	      // TLorentzVector definitions
	      TLorentzVector Lepton_1  = TLorentzVector();
	      TLorentzVector Lepton_2  = TLorentzVector();
	      TLorentzVector Lepton_3  = TLorentzVector();
	      TLorentzVector      MeT  = TLorentzVector();
	      
	      
	      Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
	      Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
	      Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index),lep_E->at(goodlep3_index));
	      MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et);
	      TLorentzVector old_nu(met_et*cos(MeT.Phi()),met_et*sin(MeT.Phi()),0.0,met_et);

  
	      float InvMass12       = (Lepton_1+Lepton_2).Mag()/1000.;
	      float InvMass13       = (Lepton_1+Lepton_3).Mag()/1000.;
	      float InvMass23       = (Lepton_2+Lepton_3).Mag()/1000.;
	      
	      float delta_lep12_Z =  -1;
	      float delta_lep13_Z =  -1;
	      float delta_lep23_Z =  -1;


	            // calculations are done step-by-step
	      
	      // charge of leptons coming from Z is different
	      if ( lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index)  < 0 ) 
		{
		  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep2_index) )
		    {
		      delta_lep12_Z =  TMath::Abs(InvMass12 - 91.18);
		    }
		}
	      
	      if ( lep_charge->at(goodlep1_index) * lep_charge->at(goodlep3_index)  < 0 ) 
		{
		  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep3_index) )
		    {
		      delta_lep13_Z =  TMath::Abs(InvMass13 - 91.18);
		    }
		}
	      
	      if ( lep_charge->at(goodlep2_index) * lep_charge->at(goodlep3_index)  < 0 ) 
		{
		  if ( lep_type->at(goodlep2_index) == lep_type->at(goodlep3_index) )
		    {
		      delta_lep23_Z =  TMath::Abs(InvMass23 - 91.18);
		    }
		}

	      

	        
	      // define candidates 
	      int wcand = 0;
	      float tmp = 0;
	      
	      
	      // begin all permutations 
	      // in the case we have eee or mumumu
	      ////		      
	      if( ( delta_lep12_Z >0 && delta_lep23_Z >0)  && delta_lep13_Z < 0 && (delta_lep12_Z < delta_lep23_Z) ) {
		tmp = delta_lep12_Z; 
		wcand = 3;
	      } 
	      
	      if( ( delta_lep12_Z >0 && delta_lep23_Z >0)   && delta_lep13_Z < 0 && (delta_lep12_Z > delta_lep23_Z) ) {
		tmp = delta_lep23_Z; 
		wcand = 1;
	      }
	      
	      if( ( delta_lep12_Z >0 && delta_lep13_Z >0 )  && delta_lep23_Z < 0 && (delta_lep12_Z < delta_lep13_Z) ) {
		tmp = delta_lep12_Z;
		wcand = 3;
	      }

	       
	      if( ( delta_lep12_Z >0 && delta_lep13_Z >0 ) && delta_lep23_Z < 0 && (delta_lep12_Z > delta_lep13_Z) ) {
		tmp = delta_lep13_Z;
		wcand = 2;
	      }
	      
	      if( ( delta_lep13_Z >0 && delta_lep23_Z >0)  && delta_lep12_Z < 0 && (delta_lep13_Z < delta_lep23_Z) ) {
		tmp = delta_lep13_Z;
		wcand = 2;
	      }
	      
	      if( ( delta_lep13_Z >0 && delta_lep23_Z >0)  && delta_lep12_Z < 0 && (delta_lep13_Z > delta_lep23_Z) ) {
		tmp = delta_lep23_Z;
		wcand = 1;
	      }


	          
	      
	      // in the case we have eemu or mumue
	      	      
	      if ( delta_lep12_Z < 0 && delta_lep23_Z  < 0 && delta_lep13_Z > 0) {
		tmp = delta_lep13_Z; 
		wcand = 2;
	      }
	      
	      if ( delta_lep12_Z < 0 && delta_lep13_Z  < 0 && delta_lep23_Z > 0) {
		tmp = delta_lep23_Z;
		wcand = 1;
	      }
	      
	      if ( delta_lep23_Z < 0 && delta_lep13_Z  < 0 && delta_lep12_Z > 0) {
		tmp = delta_lep12_Z;
		wcand = 3;
	      }

	       
	      // depending on which is the W candidate, build mtw
	      float mtw=0; float lpt; float met; float lp_x=0; float lp_y=0; float lp_z=0; float et_x=0; float et_y=0; float mu=0; float pnu_z=0;
	     
	      float mw = 80.4*1000.;
	      if(wcand==1){
		mtw = sqrt(2*Lepton_1.Pt()*MeT.Et()*(1-cos(Lepton_1.DeltaPhi(MeT))))/1000.;
		lpt = Lepton_1.Pt();
		lp_x = Lepton_1.Pt()*TMath::Cos(Lepton_1.Phi());
                lp_y = Lepton_1.Pt()*TMath::Sin(Lepton_1.Phi());
                lp_z = Lepton_1.Pt()*TMath::SinH(Lepton_1.Eta());
             
		met = MeT.Pt();
                et_x = MeT.Pt()*TMath::Cos(MeT.Phi());
                et_y = MeT.Pt()*TMath::Sin(MeT.Phi());

                mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
		float discri = (mu*mu*lp_z*lp_z - pow(Lepton_1.E()*MeT.Pt(),2) + mu*mu)/(pow(Lepton_1.Pt(),2));
                
		if(discri>=0.){
                pnu_z = 0.0;
                float check = mu*lp_z;
                if(check>0.0){
                pnu_z = mu*lp_z - TMath::Sqrt(discri) ;}
                else{
                pnu_z = mu*lp_z + TMath::Sqrt(discri) ;}
                }
                else{
                float scaling_factor = mw*mw*0.5/(met*lpt - lp_x*et_x - lp_y*et_y);
                met *= scaling_factor;
                et_x *= scaling_factor;
                et_y *= scaling_factor;
                mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
                pnu_z = mu*lp_z/(lpt*lpt);}
	      }
	      
	      if(wcand==2){
		mtw = sqrt(2*Lepton_2.Pt()*MeT.Et()*(1-cos(Lepton_2.DeltaPhi(MeT))))/1000.;
		lpt = Lepton_2.Pt();
	        lp_x = Lepton_2.Pt()*TMath::Cos(Lepton_2.Phi());
                lp_y = Lepton_2.Pt()*TMath::Sin(Lepton_2.Phi());
                lp_z = Lepton_2.Pt()*TMath::SinH(Lepton_2.Eta());
                
		met = MeT.Pt();
                et_x = MeT.Pt()*TMath::Cos(MeT.Phi());
                et_y = MeT.Pt()*TMath::Sin(MeT.Phi());

                mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
		float discri = (mu*mu*lp_z*lp_z - pow(Lepton_2.E()*MeT.Pt(),2) + mu*mu)/(pow(Lepton_2.Pt(),2));
                if(discri>=0.){
		pnu_z = 0.0;
		float check = mu*lp_z;
		if(check>0.0){
		pnu_z = mu*lp_z - TMath::Sqrt(discri) ;}
		else{
		pnu_z = mu*lp_z + TMath::Sqrt(discri) ;}
		}
		else{
		float scaling_factor = mw*mw*0.5/(met*lpt - lp_x*et_x - lp_y*et_y);
		met *= scaling_factor;
		et_x *= scaling_factor;
		et_y *= scaling_factor;
		mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
		pnu_z = mu*lp_z/(lpt*lpt);}
                }
	      
	      if(wcand==3){
		mtw = sqrt(2*Lepton_3.Pt()*MeT.Et()*(1-cos(Lepton_3.DeltaPhi(MeT))))/1000.;
		lpt = Lepton_3.Pt();
	        lp_x = Lepton_3.Pt()*TMath::Cos(Lepton_3.Phi());
                lp_y = Lepton_3.Pt()*TMath::Sin(Lepton_3.Phi());
                lp_z = Lepton_3.Pt()*TMath::SinH(Lepton_3.Eta());
  
		met = MeT.Pt();
                et_x = MeT.Pt()*TMath::Cos(MeT.Phi());
                et_y = MeT.Pt()*TMath::Sin(MeT.Phi());

                mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
               
		float discri = (mu*mu*lp_z*lp_z - pow(Lepton_3.E()*MeT.Pt(),2) + mu*mu)/(pow(Lepton_3.Pt(),2));
                if(discri>=0.){
                pnu_z = 0.0;
                float check = mu*lp_z;
                if(check>0.0){
                pnu_z = mu*lp_z - TMath::Sqrt(discri) ;}
                else{
                pnu_z = mu*lp_z + TMath::Sqrt(discri) ;}
                }
                else{
                float scaling_factor = mw*mw*0.5/(met*lpt - lp_x*et_x - lp_y*et_y);
                met *= scaling_factor;
                et_x *= scaling_factor;
                et_y *= scaling_factor;
                mu = (mw*mw)/2. + lp_x*et_x + lp_y*et_y;
                pnu_z = mu*lp_z/(lpt*lpt);}
	      }
	      
	       float pnu = pnu_z/1000.;
	       old_nu.SetPz(pnu_z);
	       //h_pnuz->Fill(pnu,weight);
	      
	      // depending on which is the W candidate, build mll and ptll, pT of lepton from W > 20 GeV     
	      float mll=0; float ptll=0; float lepW=0; float inv_mw;
	      if(wcand==1) {mll = InvMass23; ptll = (Lepton_2+Lepton_3).Pt()/1000.; lepW = Lepton_1.Pt()/1000.; pnu = (Lepton_1+old_nu).M()/1000.;}
	      if(wcand==2) {mll = InvMass13; ptll = (Lepton_1+Lepton_3).Pt()/1000.; lepW = Lepton_2.Pt()/1000.;pnu = (Lepton_2+old_nu).M()/1000.;}
	      if(wcand==3) {mll = InvMass12; ptll = (Lepton_1+Lepton_2).Pt()/1000.; lepW = Lepton_3.Pt()/1000.;pnu = (Lepton_2+old_nu).M()/1000.;}
	      

	        // cut: m_ll - m(Z) < 10
	      if(tmp < 10.)
		{
		  // mtw > 30 GeV, at least one lepton with pT > 25 GeV 
		  if( mtw > 30. && met_et > 30000. && (lep_pt->at(goodlep1_index)/1000. > 25 || lep_pt->at(goodlep2_index)/1000. > 25 || lep_pt->at(goodlep3_index)/1000. > 25) )
		  {
		    h_mtw->Fill(mtw,weight);
		    h_mll->Fill(mll,weight);
		    h_pnuz->Fill(pnu,weight);
		   
		  }
		}

	    }

       }

   }

	    
	      
	      
	
     
   //  std::cout<<"Number of Z-->mumu + 1mu events: "<<zmm_count<<std::endl;
    std::cout<<"Number of initial events: "<<count_init<<std::endl;
    std::cout<<"Number of events that passed trigger: "<<count_trig<<std::endl;
    std::cout<<"Number of nlep = 3 events: "<< count_trilep<<std::endl;
    std::cout<<"Number of nbjet = 2 ebents: "<<count_bjet<<std::endl;
    std::cout<<"Number of met passed events: "<<count_met<<std::endl;
    std::cout<<"Number of mll passed events: "<<count_mll<<std::endl;
    std::cout<<"Number of mbb passed events: "<<count_mbb<<std::endl;
    std::cout<<"Number of mllbb passed events: "<<count_mllbb<<std::endl;   
  
   
   outf.Write();
   //myfile.close();
   
}

mini::mini(TTree *tree) : fChain(0) 
 {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  if (tree == 0) {
 
      TChain* tchain = new TChain("mini");
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/Data/data_*.2lep.root");
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363356.ZqqZll.2lep.root"); // ZZ  project
      //    tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363358.WqqZll.2lep.root"); // WZ  project
      //      tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_410000.ttbar_lep.2lep.root"); // ttbar  project
      //      tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_364122.Zee_PTV140_280_BFilter.2lep.root"); // Z+jet  project
      
           tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363491.lllv.2lep.root"); // WZ  MC_actual
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363490.llll.2lep.root");  // ZZ MC
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_3641*.Zmumu_PTV0_70_CVetoBVeto.2lep.root");  //  Z + jet MC  
      //   tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392217.C1N2_WZ_400p0_0p0_3L_2L7.2lep.root"); // M1 
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
