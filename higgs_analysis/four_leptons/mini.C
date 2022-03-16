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
   TH1D* h_mllll= new TH1D("mllll","Invariant mass of 4 lep; M_{llll} (GeV/c^{2}); Events per bin",80,20.,180.);
   TH1D* h_pt= new TH1D("pt","Transverse ,momentum of 4 lep; PT_{llll} (GeV/c^{2}); Events per bin",80,20.,180.);
   TH1D* h_eta= new TH1D("eta","Rapidity of 4 lep; #eta_{llll}; Events per bin",80,-5.,5.);
  
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

  	  // Preselection of good leptons
	  int goodlep_index[lep_n];
	  int goodlep_n = 0;
	  int lep_index =0;
	  
	  for(unsigned int i=0; i<lep_n; i++)
	    {
              TLorentzVector leptemp;  leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);
	      
              // loosely isolated and very soft 
	      if( lep_pt->at(i) > 5000. && TMath::Abs(lep_eta->at(i)) < 2.5 && ( (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.3) && ( (lep_etcone20->at(i) / lep_pt->at(i)) < 0.3 ) ) {
		// electron
		if ( lep_type->at(i) == 11 && lep_pt->at(i) > 7000. && TMath::Abs(lep_eta->at(i)) <2.47 ) {
		  if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
		    goodlep_n = goodlep_n + 1;
		    goodlep_index[lep_index] = i;
		    lep_index++;
		  }
		}
		//muon
		if ( lep_type->at(i) == 13) {
		  if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 3 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
		    goodlep_n = goodlep_n + 1;
		    goodlep_index[lep_index] = i;
		    lep_index++;
		  }
		}
	      }
	    }

	  	  
	  //Exactly four good leptons
	  if(goodlep_n == 4 )
	    {
	      
	      int goodlep1_index = goodlep_index[0];
	      int goodlep2_index = goodlep_index[1];
	      int goodlep3_index = goodlep_index[2];
	      int goodlep4_index = goodlep_index[3];
	      
	      //first lepton pT > 25 GeV, second > 15 GeV and third > 10 GeV		      
	      if (lep_pt->at(goodlep1_index) > 25000. && lep_pt->at(goodlep2_index) > 15000. && lep_pt->at(goodlep3_index) > 10000. ) 
		{ 		     
		  
		  
		  // TLorentzVector definitions
		  TLorentzVector Lepton_1  = TLorentzVector();
		  TLorentzVector Lepton_2  = TLorentzVector();
		  TLorentzVector Lepton_3  = TLorentzVector();
		  TLorentzVector Lepton_4  = TLorentzVector();

		    Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
		  Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
		  Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index),lep_E->at(goodlep3_index));
		  Lepton_4.SetPtEtaPhiE(lep_pt->at(goodlep4_index), lep_eta->at(goodlep4_index), lep_phi->at(goodlep4_index),lep_E->at(goodlep4_index));
		  
		  
		  // minimisation of difference from the Z mass
		  float delta_Z1=0; 
		  float delta_Z2=0; 
		  float InvMassZ1=0; 
		  float InvMassZ2=0;
		  float delta_Z1_1=0; float delta_Z1_2=0; float delta_Z1_3=0;
		  float delta_Z2_1=0; float delta_Z2_2=0; float delta_Z2_3=0;
		  float InvMassZ1_1=0; float InvMassZ1_2=0; float InvMassZ1_3=0;
		  float InvMassZ2_1=0; float InvMassZ2_2=0; float InvMassZ2_3=0;
		  float sum_ZZ1=0; float sum_ZZ2=0; float sum_ZZ3=0;
		  
		  // final values
		  float InvMassZ1_min=0; float InvMassZ2_min=0; float sum_ZZ_fin=0;
		  
		    
		  float sum_charges = lep_charge->at(goodlep1_index) + lep_charge->at(goodlep2_index) + lep_charge->at(goodlep3_index) + lep_charge->at(goodlep4_index);			 
		  
		  // step-by-step
		  // opposite charge leptons
		  if ( sum_charges == 0  ) 
		    {
		      
		      
		      int sum_types  = lep_type->at(goodlep1_index) + lep_type->at(goodlep2_index) + lep_type->at(goodlep3_index) + lep_type->at(goodlep4_index) ;
		      
		      // type e=11, mu=13
		      // begin case e+e-e+e- or mu+mu-mu+mu-
		      if ( sum_types == 44 || sum_types == 52  )
			{
			  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep2_index) && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index)) < 0 )  )
			    {
		              InvMassZ1_1=(Lepton_1+Lepton_2).Mag()/1000.;
			      InvMassZ2_1=(Lepton_3+Lepton_4).Mag()/1000.;
			      delta_Z1_1 =  TMath::Abs(InvMassZ1_1 - 91.18); 
			      delta_Z2_1 =  TMath::Abs(InvMassZ2_1 - 91.18);
			    }
			    if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep3_index)  && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep3_index)) < 0 ) )
			    {
			      InvMassZ1_2=(Lepton_1+Lepton_3).Mag()/1000.;
			      InvMassZ2_2=(Lepton_2+Lepton_4).Mag()/1000.;
			      delta_Z1_2 =  TMath::Abs(InvMassZ1_2 - 91.18); 
			      delta_Z2_2 =  TMath::Abs(InvMassZ2_2 - 91.18);
			    }
			  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep4_index)  && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep4_index)) < 0 ) )
			    {
			      InvMassZ1_3=(Lepton_1+Lepton_4).Mag()/1000.;
			      InvMassZ2_3=(Lepton_2+Lepton_3).Mag()/1000.;
			      delta_Z1_3 =  TMath::Abs(InvMassZ1_3 - 91.18); 
			      delta_Z2_3 =  TMath::Abs(InvMassZ2_3 - 91.18);
			    }

			    if(delta_Z1_1 < delta_Z2_1) { InvMassZ1_min = InvMassZ1_1; InvMassZ2_min = InvMassZ2_1;}
                          if(delta_Z2_1 < delta_Z1_1) { InvMassZ1_min = InvMassZ2_1; InvMassZ2_min = InvMassZ1_1;}

			  if(delta_Z1_2 < delta_Z2_2) { InvMassZ1_min = InvMassZ1_2; InvMassZ2_min = InvMassZ2_2;}
                          if(delta_Z2_2 < delta_Z1_2) { InvMassZ1_min = InvMassZ2_2; InvMassZ2_min = InvMassZ1_2;}

			  if(delta_Z1_3 < delta_Z2_3) { InvMassZ1_min = InvMassZ1_3; InvMassZ2_min = InvMassZ2_3;}
                          if(delta_Z2_3 < delta_Z1_3) { InvMassZ1_min = InvMassZ2_3; InvMassZ2_min = InvMassZ1_3;}

			} // cases of eeee or mumumumu
		      
		      ////////////////////////////////////

		       // case eemumu 
		      if ( sum_types == 48 )
			{
			  
			  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep2_index)  && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_2).Mag()/1000.;
			      InvMassZ2=(Lepton_3+Lepton_4).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }
			  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep3_index)  && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep3_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_3).Mag()/1000.;
			      InvMassZ2=(Lepton_2+Lepton_4).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }
			  if ( lep_type->at(goodlep1_index) == lep_type->at(goodlep4_index)  && ( (lep_charge->at(goodlep1_index) * lep_charge->at(goodlep4_index)) < 0 ) )
			    {
			      InvMassZ1=(Lepton_1+Lepton_4).Mag()/1000.;
			      InvMassZ2=(Lepton_2+Lepton_3).Mag()/1000.;
			      delta_Z1 =  TMath::Abs(InvMassZ1 - 91.18); 
			      delta_Z2 =  TMath::Abs(InvMassZ2 - 91.18);
			    }

			    
			  if(delta_Z1 < delta_Z2) { InvMassZ1_min = InvMassZ1; InvMassZ2_min = InvMassZ2;}
			  if(delta_Z2 < delta_Z1) { InvMassZ1_min = InvMassZ2; InvMassZ2_min = InvMassZ1;}
			} // eemumu overe
		      
		      if ( (sum_types == 44 || sum_types == 52 || sum_types == 48) )
			{
			  
			  TLorentzVector FourLepSystem = TLorentzVector();
			  FourLepSystem = Lepton_1 + Lepton_2 + Lepton_3 + Lepton_4;
			  float FourLepSystem_M = FourLepSystem.Mag()/1000.;
			  float FourLepSystem_pt = FourLepSystem.Pt()/1000.;
			  float FourLepSystem_y = FourLepSystem.Rapidity();
			 
                     //Preselection of good jets
			  int goodjet_n = 0;
			  int goodjet_index = 0;
				  
			  if (jet_n > 0)
			    {
			      for(unsigned int i=0; i<jet_n; i++)
				{
				  if(jet_pt->at(i) > 30000. && TMath::Abs(jet_eta->at(i)) < 4.4)
				    {
				      goodjet_n++;
				      goodjet_index = i;
				    }
				}
			    }

			    //Start to fill histograms
			  if (goodjet_n == 0){
			    h_mllll->Fill(FourLepSystem_M,weight);
			    h_pt->Fill(FourLepSystem_pt,weight);
			    h_eta->Fill(FourLepSystem_y,weight);}
    	  


			}
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
      tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_345060.ggH125_ZZ4lep.2lep.root"); // Higgs ggh->ZZ->4l
      //    tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363491.lllv.2lep.root"); // WZ  MC_actual
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
