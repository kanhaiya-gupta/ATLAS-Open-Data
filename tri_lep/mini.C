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
 //  nentries = 10000;  // number of events

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
   TH1D* h_mwz= new TH1D("wzmass","Invariant mass of WZ boson; M_{WZ} (GeV/c^{2}); Events per bin",100,20.,3000.);
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
              
	      float Z_mass = 91.18;
	      unsigned int lW_index=-1;
	      unsigned int lZ1_index=-1;
	      unsigned int lZ2_index=-1;

	      float delta_lep12_Z =  -1;
              float delta_lep13_Z =  -1;
              float delta_lep23_Z =  -1;

	      if(lep_charge->at(goodlep1_index) + lep_charge->at(goodlep2_index) == 0.0 && lep_type->at(goodlep1_index) == lep_type->at(goodlep2_index)){
	      Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
              Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
	      delta_lep12_Z =  TMath::Abs((Lepton_1+Lepton_2).Mag()/1000. - Z_mass);}

	      if(lep_charge->at(goodlep1_index) + lep_charge->at(goodlep3_index) == 0.0 && lep_type->at(goodlep1_index) == lep_type->at(goodlep3_index)){
              Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
              Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index),lep_E->at(goodlep3_index));
              delta_lep13_Z =  TMath::Abs((Lepton_1+Lepton_3).Mag()/1000. - Z_mass);}
	      
	      if(lep_charge->at(goodlep2_index) + lep_charge->at(goodlep3_index) == 0.0 && lep_type->at(goodlep2_index) == lep_type->at(goodlep3_index)){
              Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
              Lepton_3.SetPtEtaPhiE(lep_pt->at(goodlep3_index), lep_eta->at(goodlep3_index), lep_phi->at(goodlep3_index),lep_E->at(goodlep3_index));
              delta_lep23_Z =  TMath::Abs((Lepton_2+Lepton_3).Mag()/1000. - Z_mass);}
	    

	       //define Z pair
	     float temp_diff_mass_min=-1.; 
		if( delta_lep12_Z!=Z_mass && delta_lep12_Z!=-1. ){
			temp_diff_mass_min=delta_lep12_Z;
			lZ1_index=0;
			lZ2_index=1;
			lW_index=2;
		}
		if( delta_lep13_Z!=Z_mass && delta_lep13_Z!=-1. && (temp_diff_mass_min==-1. || temp_diff_mass_min > delta_lep13_Z)){
			temp_diff_mass_min =delta_lep13_Z;
			lZ1_index=0;
			lZ2_index=2;
			lW_index=1;
		}
		if( delta_lep23_Z!=Z_mass && delta_lep23_Z!=-1. && (temp_diff_mass_min==-1. || temp_diff_mass_min > delta_lep23_Z)){
			temp_diff_mass_min =delta_lep23_Z;
			lZ1_index=1;
			lZ2_index=2;
			lW_index=0;
		}
             
		if(temp_diff_mass_min==-1.){
		continue;}

	     // TLorentzVector definitions
              TLorentzVector Lepton_W  = TLorentzVector();
              TLorentzVector Lepton_Z1  = TLorentzVector();
              TLorentzVector Lepton_Z2  = TLorentzVector();
	      TLorentzVector NeutrinoVec  = TLorentzVector();
	      TLorentzVector W  = TLorentzVector();
	      TLorentzVector Z  = TLorentzVector();
             
	      MeT.SetPtEtaPhiE(met_et, 0, met_phi , met_et); 
	      Lepton_W.SetPtEtaPhiE(lep_pt->at(lW_index), lep_eta->at(lW_index), lep_phi->at(lW_index),lep_E->at(lW_index));
	      Lepton_Z1.SetPtEtaPhiE(lep_pt->at(lZ1_index), lep_eta->at(lZ1_index), lep_phi->at(lZ1_index),lep_E->at(lZ1_index));
	      Lepton_Z2.SetPtEtaPhiE(lep_pt->at(lZ2_index), lep_eta->at(lZ2_index), lep_phi->at(lZ2_index),lep_E->at(lZ2_index));
            

              float PxMiss = MeT.Px();
              float PyMiss = MeT.Py();
              float EtMiss = MeT.Pt();
              float Wmass = 80385;
   
	      float M = Wmass * Wmass * 0.5 + Lepton_W.Px() * PxMiss + Lepton_W.Py() * PyMiss;
              double discriminant = Lepton_W.E() * Lepton_W.E() * (M * M - (EtMiss * EtMiss * Lepton_W.Perp2()));

	      	//	  real solutions (two)
              if (discriminant >= 0. ) {
              double PzMiss = 0;
              double check = M * Lepton_W.Pz();

		// We want to find the smallest solution in magnitude
              if (check > 0.0) {
                 PzMiss = (check - sqrt(discriminant)) / Lepton_W.Perp2();
	      }
	      else {
                 PzMiss = (check + sqrt(discriminant)) / Lepton_W.Perp2();
	      }

              NeutrinoVec.SetPxPyPzE(PxMiss, PyMiss, PzMiss, sqrt(EtMiss * EtMiss + PzMiss * PzMiss));
           }
	  
	     //	  imaginary solutions
	    else {

                double scaling_factor = Wmass * Wmass * 0.5 / (EtMiss * Lepton_W.Perp() - Lepton_W.Px() * PxMiss - Lepton_W.Py() * PyMiss);
                EtMiss *= scaling_factor;
                PxMiss *= scaling_factor;
                PyMiss *= scaling_factor;
                M = Wmass * Wmass * 0.5 + Lepton_W.Px() * PxMiss + Lepton_W.Py() * PyMiss;
                double PzMiss = (Lepton_W.Pz() * M) / Lepton_W.Perp2();
                NeutrinoVec.SetPxPyPzE(PxMiss, PyMiss, PzMiss, sqrt(EtMiss * EtMiss + PzMiss * PzMiss));
              }
          

	    W.SetPxPyPzE(NeutrinoVec.Px()+Lepton_W.Px(),NeutrinoVec.Py()+Lepton_W.Py(),NeutrinoVec.Pz()+Lepton_W.Pz(),NeutrinoVec.E()+Lepton_W.E());
            Z.SetPxPyPzE(Lepton_Z1.Px()+Lepton_Z2.Px(),Lepton_Z1.Py()+Lepton_Z2.Py(),Lepton_Z1.Pz()+Lepton_Z2.Pz(),Lepton_Z1.E()+Lepton_Z2.E());
      
	  //  std::cout<<"test10 "<<std::endl;

	    double mtw = sqrt(2*Lepton_W.Pt()*MeT.Et()*(1-cos(Lepton_W.DeltaPhi(MeT))))/1000.;
	    double mw = W.M()/1000.;
            double mll = Z.M()/1000.;
	    double mwz = (W + Z).M()/1000.;
            h_wmass->Fill(mw,weight); 
            h_mtw->Fill(mtw,weight);
            h_mll->Fill(mll,weight);
	    h_mwz->Fill(mwz,weight);

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
