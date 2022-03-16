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
   int count_lep = 0;
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
   TH1D* h_leppt= new TH1D("leppt","Leading P_{T} lepton; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   TH1D* h_subleppt= new TH1D("subleppt","Sub-leading P_{T} lepton; Sub lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   TH1D* h_jetpt= new TH1D("jetpt","Leading P_{T} b-jet; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   TH1D* h_subjetpt= new TH1D("subjetpt","Sub-leading P_{T} b-jet; Lead P_{T} (GeV/c); Events per bin", 200, 0., 500.);
   TH1D* h_dilep= new TH1D("dilep","Di-lepton invariant mass; M_{ll} (GeV/c^{2}); Events per bin",200,0.,800.);
   TH1D* h_dibjet= new TH1D("dibjet","Di b-jet invariant mass; M_{bb} (GeV/c^{2}); Events per bin",200,0.,800.);
   TH1D* h_lepjet= new TH1D("lepjet","Lepton b-jet invariant mass; M_{llbb} (GeV/c^{2}); Events per bin",200,0.,2000.);

   TH2D* h_lepjet2d = new TH2D("lepjet2d","2D Histogram plot ;M_{ll}/92 (GeV/c^{2});M_{bb}/125 (GeV/c^{2})",100,0.6,1.4,100, 0.2, 2.);

   //   TH1D* h_leadlepmass= new TH1D("leadlepmass","Lead lep Pt invariant mass; M_{ll}^{lead}; Events per bin",200,0.,500.);
   
   // TH1D* h_mupt= new TH1D("mupt","Muon Transverse Momentum; P_{T} (GeV/c); Events per bin",200,0.,100.);
   // TH1D* h_mdimu= new TH1D("mass_dimu","Di-muon invariant mass; M_{#mu^{+}#mu^{-}}(GeV/c^{2});Events per bin",400,0.,200.);
   // TH1D* h_mtrimu= new TH1D("mass_trimu","Tri-muon invariant mass; M (GeV/c^{2});Events per bin",200,0.,2000.);
  
  
   
   for  (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      count_init = count_init + 1;

      
    //Scale factors
    // luminosity = N/XSection = nentries*mc/XSection;
    // To make the luminosity same, sacle = (lumi of data)/(lumi of MC);
      
     float scaleFactor = 6.181*scaleFactor_ELE*scaleFactor_MUON*scaleFactor_LepTRIGGER*scaleFactor_PILEUP;
      
      //MC weight
      float  m_mcWeight = mcWeight;

      //Total weight
      //   float weight = scaleFactor*m_mcWeight*XSection/1000.;

        float weight = 1.;    // data


      
//Preselection cut for electron/muon trigger
	if(trigE || trigM){

  count_trig = count_trig + 1;
	  
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
     count_lep = count_lep + 1;   
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
	   count_bjet = count_bjet + 1;
	   
	      if (met_et<30000.)
            {
	      count_met = count_met + 1;

           int goodbjet1_index = goodbjet_index[0];
	   int goodbjet2_index = goodbjet_index[1];
                         
           // TLorentzVector definitions
	   TLorentzVector bjet_1  = TLorentzVector();
	   TLorentzVector bjet_2  = TLorentzVector();

        bjet_1.SetPtEtaPhiE(jet_pt->at(goodbjet1_index), jet_eta->at(goodbjet1_index), jet_phi->at(goodbjet1_index),jet_E->at(goodbjet1_index));
	bjet_2.SetPtEtaPhiE(jet_pt->at(goodbjet2_index), jet_eta->at(goodbjet2_index), jet_phi->at(goodbjet2_index),jet_E->at(goodbjet2_index));
				  
       float Mjjmax= ( bjet_1 + bjet_2 ).M()/1000.; // first indices
       float lepjet_mass =  ( Lepton_1 + Lepton_2 + bjet_1 + bjet_2 ).M()/1000.;

       h_dilep->Fill(InvMass_Leptons,weight);
       h_dibjet->Fill(Mjjmax,weight);
       h_lepjet->Fill(lepjet_mass,weight);

       h_lepjet2d->Fill(InvMass_Leptons/92.,Mjjmax/125.);

       if (InvMass_Leptons>= 80. && InvMass_Leptons <=100.)
	 { count_mll = count_mll + 1;
	   if (Mjjmax>=100. && Mjjmax <=140.)
	     { count_mbb = count_mbb + 1;
	   if (lepjet_mass>=600. && lepjet_mass<=1000)
	     { count_mllbb = count_mllbb + 1;}}}

       // How can I add this to the histogram
       //   TF2 *f2 = new TF2("f2","((x-1)/0.1)^2 + ((y-1)/0.05)^2 = 1.",0,5,0,5);
       // f2->Draw("same");


       //finding the leading and sub-leading pt leptons
			  
       double names_of_leadlep_variable[]={Lepton_1.Pt()/1000., Lepton_1.Eta(), Lepton_1.E()/1000., Lepton_1.Phi(), (double)lep_charge->at(goodlep1_index), (double)lep_type->at(goodlep1_index)};
			  
       double names_of_subleadlep_variable[]={Lepton_2.Pt()/1000., Lepton_2.Eta(), Lepton_2.E()/1000., Lepton_2.Phi(),(double)lep_charge->at(goodlep2_index), (double)lep_type->at(goodlep2_index)};

       //Start to fill histograms : find the histogram array length
       int length_leadlep = sizeof(names_of_leadlep_variable)/sizeof(names_of_leadlep_variable[0]);
       int length_subleadlep = sizeof(names_of_subleadlep_variable)/sizeof(names_of_subleadlep_variable[0]);
			  
       //Fill histograms
       for (int i=0; i<length_leadlep; i++)
	 { if (names_of_leadlep_variable[i] > 25.)
	     {h_leppt->Fill(names_of_leadlep_variable[i],weight);
	       h_subleppt->Fill(names_of_subleadlep_variable[i],weight); }
	 }

       //finding the leading and sub-leading trnsverse momentum b-jets

       double names_of_leadjet_variable[]={bjet_1.Pt()/1000., bjet_1.Eta(), bjet_1.E()/1000., Lepton_1.Phi(), (double)jet_MV2c10->at(goodbjet1_index)};
        double names_of_subleadjet_variable[]={bjet_2.Pt()/1000., bjet_2.Eta(), bjet_2.E()/1000., Lepton_2.Phi(), (double)jet_MV2c10->at(goodbjet1_index)};

       int length_leadbjet = sizeof(names_of_leadjet_variable)/sizeof(names_of_leadlep_variable[0]);


       int length_subleadbjet = sizeof(names_of_subleadjet_variable)/sizeof(names_of_subleadlep_variable[0]);

       // Fill histograms
       for (int i =0; i<length_leadbjet; i++)
	 {if (names_of_leadjet_variable[i]> 30.)
	     { h_jetpt->Fill(names_of_leadjet_variable[i],weight);
	       h_subjetpt->Fill(names_of_subleadjet_variable[i],weight);}
	 }
       
	                 
       //      if(type_one==11) FillHistogramsGlobal(InvMass_Leptons_ee , weight, "hist_ee_mLL");
       //    if(type_one==13) FillHistogramsGlobal(InvMass_Leptons_mumu , weight, "hist_mumu_mLL");
       //    FillHistogramsGlobal(InvMass_Leptons, weight, "hist_mLL");

 

       }//jet cut 


			}
		    }
		}
	    }
	}
 
  
   }
   }
  
   //  std::cout<<"Number of Z-->mumu + 1mu events: "<<zmm_count<<std::endl;
    std::cout<<"Number of initial events: "<<count_init<<std::endl;
    std::cout<<"Number of events that passed trigger: "<<count_trig<<std::endl;
    std::cout<<"Number of nlep = 2 events: "<<count_lep<<std::endl;
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
          tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/Data/data_*.2lep.root");
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363356.ZqqZll.2lep.root"); // ZZ  project
      //    tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363358.WqqZll.2lep.root"); // WZ  project
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_410000.ttbar_lep.2lep.root"); // ttbar  project
      //      tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_364122.Zee_PTV140_280_BFilter.2lep.root"); // Z+jet  project
      
      //     tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363491.lllv.2lep.root"); // WZ  MC_actual
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
