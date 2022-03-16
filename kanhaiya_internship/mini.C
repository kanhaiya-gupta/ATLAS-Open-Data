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
  
  
   
   for  (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
     
 
//Preselection cut for electron/muon trigger
 if(trigE || trigM)
{
	  
  // Preselection of good leptons
 int goodlep_index[2];
 int goodlep_n = 0;
 int lep_index =0;
	  
  for(unsigned int i=0; i<lep_n; i++)
    {
        // temporary
        TLorentzVector leptemp;
        leptemp.SetPtEtaPhiE(lep_pt->at(i)/1000., lep_eta->at(i), lep_phi->at(i), lep_E->at(i)/1000.);  


       // Lepton is Tight
       if( lep_isTightID->at(i) )
	  {
	    // Lepton is isolated and hard pT
	 if( lep_pt->at(i) >25000. && ( (lep_ptcone30->at(i)/lep_pt->at(i)) < 0.15) && ( (lep_etcone20->at(i) / lep_pt->at(i)) < 0.15 ) )
	   {
// electron selection in fiducial region excluding candidates in the transition region between the barrel and endcap electromagnetic calorimeters
      if ( lep_type->at(i) == 11 && TMath::Abs(lep_eta->at(i)) < 2.47 && ( TMath::Abs(lep_eta->at(i)) < 1.37 || TMath::Abs(lep_eta->at(i)) > 1.52 ) ) 
      {
     if( TMath::Abs(lep_trackd0pvunbiased->at(i))/lep_tracksigd0pvunbiased->at(i) < 5 && TMath::Abs(lep_z0->at(i)*TMath::Sin(leptemp.Theta())) < 0.5) {
       	            goodlep_n = goodlep_n + 1;
		    goodlep_index[lep_index] = i;
		    lep_index++;
		         }
                      }
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

 
  
  
		  
 //Exactly two good lepton
     if(goodlep_n==2)
         {
		      
	    //Preselection of good jets
	    int goodbjet_index[2];
            int goodbjet_n = 0;
	    int bjet_index = 0;
		      
		     
		      
	  for(unsigned int i=0; i<jet_n; i++)
             {
	        if(jet_pt->at(i) > 30000. && TMath::Abs(jet_eta->at(i)) < 2.5)
		   {
			 // JVT cleaning
			bool jvt_pass=true;
		        if (jet_pt->at(i) < 60000. && TMath::Abs(jet_eta->at(i)) < 2.4 && jet_jvt->at(i) < 0.59) jvt_pass=false;
		        if (jvt_pass) 
			     {			  
				// cut on 0.8244273 is 70% WP	
				if (jet_MV2c10->at(i) >0.8244273)
				   {
				      goodbjet_n = goodbjet_n + 1;
				      goodbjet_index[bjet_index] = i;
				      bjet_index++;
				    }
				}
			    }
			}

	
	      
   int goodlep1_index = goodlep_index[0];
   int goodlep2_index = goodlep_index[1];
	      
   // TLorentzVector definitions
   TLorentzVector Lepton_1  = TLorentzVector();
   TLorentzVector Lepton_2  = TLorentzVector();
	      
  Lepton_1.SetPtEtaPhiE(lep_pt->at(goodlep1_index), lep_eta->at(goodlep1_index), lep_phi->at(goodlep1_index),lep_E->at(goodlep1_index));
  Lepton_2.SetPtEtaPhiE(lep_pt->at(goodlep2_index), lep_eta->at(goodlep2_index), lep_phi->at(goodlep2_index),lep_E->at(goodlep2_index));
	      
	      
  TLorentzVector     Lepton_12 = TLorentzVector();
  Lepton_12 = Lepton_1 + Lepton_2;
  float InvMass_Leptons = Lepton_12.Mag()/1000.;
 
   
//Leptons of opposite charge
if(lep_charge->at(goodlep1_index) * lep_charge->at(goodlep2_index)  < 0)
 {
   // Leptons of same flavour
   int type_one = lep_type->at(goodlep1_index);
   int type_two = lep_type->at(goodlep2_index);
   if(TMath::Abs(type_one) == TMath::Abs(type_two))
   {
     float InvMass_Leptons_ee = 0.; if(type_one==11) InvMass_Leptons_ee = InvMass_Leptons;
     float InvMass_Leptons_mumu = 0.; if(type_one==13) InvMass_Leptons_mumu = InvMass_Leptons;
		      
     // Invariant mass selection: m_ll - mZ < 25 GeV
     if( (TMath::Abs(InvMass_Leptons_ee - 91.18) < 25. ) || (TMath::Abs(InvMass_Leptons_mumu - 91.18) < 25. ) )
	{
			  
// By default, we are using for this analysis a MC sample known to describe poorly large jet multiplicity, thus we cut on nJets==0, lepton kinematics are well described in this phase-space 
			  //    FillHistogramsLeadJet((double)jet_n, weight, "hist_n_jets");
	  
      if(goodbjet_n==2)
         {
           int goodbjet1_index = goodbjet_index[0];
	   int goodbjet2_index = goodbjet_index[1];
                         
           // TLorentzVector definitions
	   TLorentzVector bjet_1  = TLorentzVector();
	   TLorentzVector bjet_2  = TLorentzVector();

        bjet_1.SetPtEtaPhiE(jet_pt->at(goodbjet1_index), jet_eta->at(goodbjet1_index), jet_phi->at(goodbjet1_index),jet_E->at(goodbjet1_index));
	bjet_2.SetPtEtaPhiE(jet_pt->at(goodbjet2_index), jet_eta->at(goodbjet2_index), jet_phi->at(goodbjet2_index),jet_E->at(goodbjet2_index));
				  
       float Mjjmax= ( bjet_1 + bjet_2 ).M()/1000.; // first indices  
	                 
       //      if(type_one==11) FillHistogramsGlobal(InvMass_Leptons_ee , weight, "hist_ee_mLL");
       //    if(type_one==13) FillHistogramsGlobal(InvMass_Leptons_mumu , weight, "hist_mumu_mLL");
       //    FillHistogramsGlobal(InvMass_Leptons, weight, "hist_mLL");

 

       }//jet cut 


			}
		    }
		}
	    }
	}
 
  
      
    

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

    
      

      
     
      
   
   }

  
   //  std::cout<<"Number of Z-->mumu + 1mu events: "<<zmm_count<<std::endl;
   
  
   
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
