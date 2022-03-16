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
   TH1D* h_mtrimu= new TH1D("mass_trimu","Tri-muon invariant mass; M (GeV/c^{2});Events per bin",400,0.,600.);
   TH1D* h_deltaR= new TH1D("deltaR","deltaR;#Delta R; Events per bin",400,0.,10.0);
   TH1D* h_charge= new TH1D("sum_charge","sum_charge; #Sigma_{q_{i}};Events per bin",40,-10,10);
   TH1D* h_ptcone1 = new TH1D("ptconemu1","ptconemu1", 400,0.,0.8);
   TH1D* h_ptcone2 = new TH1D("ptconemu2","ptconemu2", 400,0.,0.8);
   TH1D* h_ptcone3 = new TH1D("ptconemu3","ptconemu3", 400,0.,0.8);
   TH1D* h_met = new TH1D("met","met; Missing Transverse Energy; #slash{E_{T}} (GeV); Events per bin",100,0.,600.);
   TH1D* h_Mt = new TH1D("M_T(W)","Transverse mass ;M_{T}[W] (GeV/c^{2});Events per bin",25,0.,200.);
   TH1D* h_Mt_wz = new TH1D("M_T(WZ)","Transverse mass of WZ; M_{T}[WZ](GeV/c^{2}); Events per bin",80,40.,1000.);
   TH1D* h_M_wz = new TH1D("M_(WZ)","Invariant mass of WZ; M_{inv}[WZ](GeV/c^{2});Events per bin",80,40.,1000.);
   TH1D* h_nan = new TH1D("N_nan","Histogram plot of  the argument of cosh^-1(..) BSM3 WZ_500p;;Events per bin", 80,-10., 30.);
   TH2D* h2_plot1 = new TH2D("2d_plot1","Et_mis vs M_T(WZ);M_T(WZ);Et_mis",50,0.,800.,50, 10., 200.);
   TH2D* h2_plot2 = new TH2D("2d_plot2","",50,0.,800.,50,10.,200.);

   TH2D* h_mtw_tri = new TH2D("2d_mtw_tri","2D Plot of M_{T}(W) vs Tri-muons inv mass;M_{T}(W);M_{3#muon}(inv)",100,40.,800.,100,40.,200.);

   TH1D* h_phi = new TH1D("d_phi","Cos(#Delta #Phi) plot of #mu #nu;;Events per bin",80,-5.,5.);
   TH1D* h_eta = new TH1D("d_eta","#Delta #eta plot;;Events per bin",80,-10.,10.);
   TH1D* h_arg = new TH1D("d_arg","m_w^2/(2*Pt*Etmis);;Events per bin",80,-2.,20.);


    TH1D* h_Mt_weit = new TH1D("M_T_weit(W)","Transverse mass ;M_{T}[W] (GeV/c^{2});Events per bin",25,0.,200.);
    TH1D* h_Mt_wz_weit = new TH1D("M_T_weit(WZ)","Transverse mass of WZ; M_{T}[WZ](GeV/c^{2}); Events per bin",80,0.,1000.);

    TH1D* H_Mt_weit_mc = new TH1D("M_T_WEIT-MC(W)","TRANSVERSE MASS ;M_{T}[W] (GEV/C^{2});EVENTS PER BIN",25,0.,200.);
    TH1D* H_MT_WZ_WEIT_MC = new TH1D("M_T_WEIT_MC(WZ)","TRANSVERSE MASS OF WZ; M_{T}[WZ](GEV/C^{2}); EVENTS PER BIN",80,40.,1000.);

    TH1D* H_DELTARNU = new TH1D("DELTA_RNU","#DELTA#PHI BETWEEN #MU & #NU;#DELTA#PHI; EVENTS PER BIN",80,-10,10.);
   
   // DI-MONS HISTOGRAMS
   /*  
   TH1D* D1_MUPT= NEW TH1D("D1-MUPT","D1-MUPT",200,0.,100.);
   TH1D* D2_MUPT= NEW TH1D("D2-MUPT","D2-MUPT",200,0.,100.);
   TH1D* D_MDIMU= NEW TH1D("D-MASS_DIMU","D-MASS_DIMU",400,0.,200.);
   TH1D* D1_MT = NEW TH1D("D1-M_T","D1-M_T",400,0.,200.);
   TH1D* D2_MT = NEW TH1D("D2-M_T","D2-M_T",400,0.,200.);
   TH1D* DZ_MT = NEW TH1D("DZ-M_T","DZ-M_T",400,0.,200.);
   */ 
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      // CALCULATE M?T FOR THE FIRST MUON IN DI-LEPTON EVENTS
      /* 
      IF (LEP_N != 2) CONTINUE;
      IF ((*LEP_TYPE)[0] != 11) CONTINUE; // SELECT 3 MU
      IF ((*LEP_TYPE)[1] != 13) CONTINUE; // SELECT 3 MU
      IF ((*LEP_PT)[0]>13.E6 OR (*LEP_PT)[1]>13.E6) CONTINUE; // DISCARD MISMEASURED EVENTS
      IF ((*LEP_ISTIGHTID)[0] == FALSE) CONTINUE; 
      IF ((*LEP_ISTIGHTID)[1] == FALSE) CONTINUE;
       IF ((*LEP_PT)[0]<20000.) CONTINUE;
      IF ((*LEP_PT)[1]<20000.) CONTINUE;

      TLORENTZVECTOR DMU1,DMU2,MET1; // 4-MOMENTUM IN GEV
	    
      DMU1.SETPTETAPHIM((*LEP_PT)[0]/1000.,(*LEP_ETA)[0],(*LEP_PHI)[0],0.1057);
      DMU2.SETPTETAPHIM((*LEP_PT)[1]/1000.,(*LEP_ETA)[1],(*LEP_PHI)[1],0.0005);
      MET1.SETPTETAPHIE((MET_ET)/1000.,0.,MET_PHI,(MET_ET)/1000.);

      FLOAT DMU12 = (DMU1 + DMU2).M();
       
       D1_MUPT->FILL(DMU1.PT());
       D2_MUPT->FILL(DMU2.PT());
       D_MDIMU->FILL(DMU12);

       FLOAT MTW1 = 0;
       FLOAT MTW2 = 0;
       FLOAT MTZ = 0;

       IF (TMATH::ABS(DMU12 - 92.1) > 15)
      {  MTW1 = SQRT(2*DMU1.PT()*MET1.E()*(1-COS(DMU1.DELTAPHI(MET1))));
       MTW2 = SQRT(2*DMU2.PT()*MET1.E()*(1-COS(DMU2.DELTAPHI(MET1))));
       MTZ = SQRT(2*DMU1.PT()*DMU2.PT()*(1-COS(DMU1.DELTAPHI(DMU2))));
       
       
       D1_MT->FILL(MTW1);
       D2_MT->FILL(MTW2);
       DZ_MT->FILL(MTZ);
      }

       // DI MUON END

       */
     
      
      				 
      // SELECT E MU FINAL STATE
      IF (LEP_N != 3) CONTINUE;
      IF ((*LEP_TYPE)[0] != 13) CONTINUE; // SELECT 3 MU
      IF ((*LEP_TYPE)[1] != 13) CONTINUE; // SELECT 3 MU
      IF ((*LEP_TYPE)[2] != 13) CONTINUE; // SELECT 3 MU
      IF ((*LEP_PT)[0]>13.E6 OR (*LEP_PT)[1]>13.E6 OR (*LEP_PT)[2]>13.E6) CONTINUE; // DISCARD MISMEASURED EVENTS

      // TIGHT LEPTON SELECTION
      IF ((*LEP_ISTIGHTID)[0] == FALSE) CONTINUE; 
      IF ((*LEP_ISTIGHTID)[1] == FALSE) CONTINUE;
      IF ((*LEP_ISTIGHTID)[2] == FALSE) CONTINUE;

      // HARD LEPTONS
      IF ((*LEP_PT)[0]<20000.) CONTINUE;
      IF ((*LEP_PT)[1]<20000.) CONTINUE;
      IF ((*LEP_PT)[2]<20000.) CONTINUE;
      
      TLORENTZVECTOR PMU1,PMU2,PMU3,MET,PNU1, PNU2; // 4-MOMENTUM IN GEV
      PMU1.SETPTETAPHIM((*LEP_PT)[0]/1000.,(*LEP_ETA)[0],(*LEP_PHI)[0],0.1057);
      PMU2.SETPTETAPHIM((*LEP_PT)[1]/1000.,(*LEP_ETA)[1],(*LEP_PHI)[1],0.1057);
      PMU3.SETPTETAPHIM((*LEP_PT)[2]/1000.,(*LEP_ETA)[2],(*LEP_PHI)[2],0.1057);
      MET.SETPTETAPHIM((MET_ET)/1000.,0.,MET_PHI,0.);
      
      FLOAT M2MU = 0;
      FLOAT M3MU = (PMU1+PMU2+PMU3).M();
      FLOAT DELTA_M1M2 = PMU1.DELTAR(PMU2);
      FLOAT DELTA_M1M3 = PMU1.DELTAR(PMU3);
      FLOAT DELTA_M2M3 = PMU2.DELTAR(PMU3);
      FLOAT MAX_R = TMATH::MAX(DELTA_M1M2, DELTA_M1M3);
      FLOAT MAX_DELR = TMATH::MAX(MAX_R, DELTA_M2M3);
      INT LEP_CHAR = (*LEP_CHARGE)[0] + (*LEP_CHARGE)[1] + (*LEP_CHARGE)[2];
      FLOAT PT_CONE1 = (*LEP_PTCONE30)[0]/(*LEP_PT)[0];
      FLOAT PT_CONE2 = (*LEP_PTCONE30)[1]/(*LEP_PT)[1];
      FLOAT PT_CONE3 = (*LEP_PTCONE30)[2]/(*LEP_PT)[2];

      //SCALE FACTORS
      // LUMINOSITY = N/XSECTION = NENTRIES*MC/XSECTION;
      // TO MAKE THE LUMINOSITY SAME, SACLE = (LUMI OF DATA)/(LUMI OF MC);
      
      FLOAT LUMI_MC = 10.*XSECTION*1000./(NENTRIES*MCWEIGHT);
      FLOAT SCALEFACTOR = SCALEFACTOR_ELE*SCALEFACTOR_MUON*SCALEFACTOR_LEPTRIGGER*SCALEFACTOR_PILEUP*XSECTION/1000;
      
      //MC WEIGHT
      FLOAT_T M_MCWEIGHT = MCWEIGHT;

      //TOTAL WEIGHT
      //  FLOAT WEIGHT = 7.8*SCALEFACTOR*M_MCWEIGHT;    // XSECTION/1000.;

      FLOAT WEIGHT = 0.15*M_MCWEIGHT;   // FOR DATA //  0.0025; 

      // TRANSVERSE MSS OF W BOSON CALCULATION

      FLOAT M12 = (PMU1 + PMU2).M();
      FLOAT M13 = (PMU1 + PMU3).M();
      FLOAT M23 = (PMU2 + PMU3).M();

      FLOAT DELTA_M12 =-1.;
      FLOAT DELTA_M13 =-1.;
      FLOAT DELTA_M23 =-1.;
     

      // LEPTONS COMING FROM THE Z HAVE DIFFERENT CHARGE

      IF ( (*LEP_CHARGE)[0]* (*LEP_CHARGE)[1] < 0)
	{ DELTA_M12 = TMATH::ABS(M12-91.2);}

       IF ( (*LEP_CHARGE)[0]* (*LEP_CHARGE)[2] < 0)
	{ DELTA_M13 = TMATH::ABS(M13-91.2);}

        IF ( (*LEP_CHARGE)[1]* (*LEP_CHARGE)[2] < 0)
	{ DELTA_M23 = TMATH::ABS(M23-91.2);}

      // DEFINE CANDIDATES 
	 INT WCAND = 0;
	 FLOAT TMP = 0.;
	

	  // IN THE CASE WE HAVE EEE OR MUMUMU
	      ////		      
	      IF( ( DELTA_M12 >0 && DELTA_M23 >0)  && DELTA_M13 < 0 && (DELTA_M12 < DELTA_M23) ) {
		TMP = DELTA_M12; 
		WCAND = 3;
		ZMM_COUNT++;
	      } 
	      
	      IF( ( DELTA_M12 >0 && DELTA_M23 >0)   && DELTA_M13 < 0 && (DELTA_M12 > DELTA_M23) ) {
		TMP = DELTA_M23; 
		WCAND = 1;
		ZMM_COUNT++;
	      }
	      
	      IF( ( DELTA_M12 >0 && DELTA_M13 >0 )  && DELTA_M23 < 0 && (DELTA_M12 < DELTA_M13) ) {
		TMP = DELTA_M12;
		WCAND = 3;
		ZMM_COUNT++;
	      }
	      
	      IF( ( DELTA_M12 >0 && DELTA_M13 >0 ) && DELTA_M23 < 0 && (DELTA_M12 > DELTA_M13) ) {
		TMP = DELTA_M13;
		WCAND = 2;
		ZMM_COUNT++;
	      }
	      
	      IF( ( DELTA_M13 >0 && DELTA_M23 >0)  && DELTA_M12 < 0 && (DELTA_M13 < DELTA_M23) ) {
		TMP = DELTA_M13;
		WCAND = 2;
		ZMM_COUNT++;
	      }
	      
	      IF( ( DELTA_M13 >0 && DELTA_M23 >0)  && DELTA_M12 < 0 && (DELTA_M13 > DELTA_M23) ) {
		TMP = DELTA_M23;
		WCAND = 1;
		ZMM_COUNT++;
	      }
	 

      FLOAT MTW = 0.;
      FLOAT ETA_1 = 0.;
      FLOAT ETA_2 = 0.;
      FLOAT MWZ = 0.;
      FLOAT ARG_COSH = 0.;
      FLOAT D_PHI = 0;
      FLOAT D_ETA = 0.;
      FLOAT ARG_1 = 0.;

       IF  (WCAND == 1) 
       	{  MTW = SQRT(2*PMU1.PT()*MET.PT()*(1-COS(PMU1.DELTAPHI(MET))));
	   M2MU =  (PMU2+PMU3).M();
	  ETA_1 = PMU1.ETA() - TMATH::ACOSH((MW*MW)/(2*PMU1.PT()*MET.E()) + COS(PMU1.DELTAPHI(MET)));
	  ETA_2 = PMU1.ETA() + TMATH::ACOSH((MW*MW)/(2*PMU1.PT()*MET.E()) + COS(PMU1.DELTAPHI(MET)));
	  PNU1.SETPTETAPHIM((MET_ET)/1000.,ETA_1,MET_PHI,0.);
	  PNU2.SETPTETAPHIM((MET_ET)/1000.,ETA_2,MET_PHI,0.);
	  MWZ = (PMU1 + PMU2 +  PMU3 + PNU1).M();
          ARG_COSH =  COS(PMU1.DELTAPHI(MET)); // (MW*MW)/(2*PMU1.PT()*MET.E()) + COS(PMU1.DELTAPHI(MET));

	  D_PHI = COS(PMU1.DELTAPHI(MET));
	  D_ETA = PMU1.ETA() - MET.ETA();
	  ARG_1 = (MW*MW)/(2*PMU1.PT()*MET.E());

	   H_DELTARNU->FILL(PMU1.DELTAPHI(MET));
	  
	  //  MWZ = (PMU1 + PMU2 +  PMU3 + PNU2).M();
	   
	}

        IF  (WCAND == 2)
       	{  MTW = SQRT(2*PMU2.PT()*MET.PT()*(1-COS(PMU2.DELTAPHI(MET))));
	   M2MU =  (PMU1+PMU3).M();
	   ETA_1 = PMU2.ETA() - TMATH::ACOSH((MW*MW)/(2*PMU2.PT()*MET.E()) + COS(PMU2.DELTAPHI(MET)));
           ETA_2 = PMU2.ETA() + TMATH::ACOSH((MW*MW)/(2*PMU2.PT()*MET.E()) + COS(PMU2.DELTAPHI(MET)));
	   PNU1.SETPTETAPHIM((MET_ET)/1000.,ETA_1,MET_PHI,0.);
	   PNU2.SETPTETAPHIM((MET_ET)/1000.,ETA_2,MET_PHI,0.);
	   MWZ = (PMU1 + PMU2 +  PMU3 + PNU1).M();
	   ARG_COSH =  COS(PMU2.DELTAPHI(MET)); // (MW*MW)/(2*PMU1.PT()*MET.E()) + COS(PMU2.DELTAPHI(MET));

	   D_PHI = COS(PMU2.DELTAPHI(MET));
	   D_ETA = PMU2.ETA() - MET.ETA();
	   ARG_1 = (MW*MW)/(2*PMU2.PT()*MET.E());
	   //   MWZ = (PMU1 + PMU2 +  PMU3 + PNU2).M();
	   H_DELTARNU->FILL(PMU2.DELTAPHI(MET));
	}
       
      IF  (WCAND == 3)
	 {  MTW = SQRT(2*PMU3.PT()*MET.PT()*(1-COS(PMU3.DELTAPHI(MET))));
	    M2MU =  (PMU1+PMU2).M();
            ETA_1 = PMU3.ETA() - TMATH::ACOSH((MW*MW)/(2*PMU3.PT()*MET.E()) + COS(PMU3.DELTAPHI(MET)));
	    ETA_2 = PMU3.ETA() + TMATH::ACOSH((MW*MW)/(2*PMU3.PT()*MET.E()) + COS(PMU3.DELTAPHI(MET)));
	    PNU1.SETPTETAPHIM((MET_ET)/1000.,ETA_1,MET_PHI,0.);
	    PNU2.SETPTETAPHIM((MET_ET)/1000.,ETA_2,MET_PHI,0.);
	    MWZ = (PMU1 + PMU2 +  PMU3 + PNU1).M();
	    ARG_COSH = COS(PMU1.DELTAPHI(MET)); //  (MW*MW)/(2*PMU1.PT()*MET.E()) + COS(PMU1.DELTAPHI(MET));

	    D_PHI = COS(PMU1.DELTAPHI(MET));
	    D_ETA = PMU3.ETA() - MET.ETA();
	    ARG_1 = (MW*MW)/(2*PMU3.PT()*MET.E());
	    //  MWZ = (PMU1 + PMU2 +  PMU3 + PNU2).M();

	    H_DELTARNU->FILL(PMU3.DELTAPHI(MET));
	 }


      // CUT: M_LL - M(Z) < 10
	      IF(TMP < 10.)
		{
		  // MTW > 30 GEV, AT LEAST ONE LEPTON WITH PT > 25 GEV 
		  IF( MTW > 30. && MET.PT() > 30. && (PMU1.PT() > 25. || PMU2.PT() > 25. || PMU3.PT() > 25.) )
		    {
		      H_MT->FILL(MTW);
		      H_MT_WEIT->FILL(MTW,WEIGHT);
		       H_MTW_TRI->FILL(MTW, M3MU);
		       H_MT_WEIT_MC->FILL(MTW,WEIGHT);
		      }}
	     
	     
	      H_M_WZ->FILL(MWZ,WEIGHT);
	      H_MUPT->FILL(PMU1.PT(),WEIGHT);
	      H_MUPT->FILL(PMU2.PT(),WEIGHT);
	      H_MUPT->FILL(PMU3.PT(),WEIGHT);
	      H_MDIMU->FILL(M2MU,WEIGHT);
      // IF (DELTA_M1M2 < 0.4 && DELTA_M1M3 < 0.4 && DELTA_M2M3 < 0.4)
      IF (MAX_DELR > 0.4)
	{ H_MTRIMU->FILL(M3MU,WEIGHT);}
      H_DELTAR->FILL(MAX_DELR,WEIGHT);
      H_CHARGE->FILL(LEP_CHAR,WEIGHT);
      H_PTCONE1->FILL(PT_CONE1,WEIGHT);
      H_PTCONE2->FILL(PT_CONE2,WEIGHT);
      H_PTCONE3->FILL(PT_CONE3,WEIGHT);
      H_MET->FILL(MET.PT(),WEIGHT);
      H_NAN->FILL(ARG_COSH,WEIGHT);

      H_PHI->FILL(D_PHI,WEIGHT);
      H_ETA->FILL(D_ETA,WEIGHT);
      H_ARG->FILL(ARG_1,WEIGHT);

     
      // CALCULATION FOR THE TRANSVERSE MASS OF WZ BOSON PAIR 

      FLOAT P_X1 = PMU1.PT()*TMATH::COS(PMU1.PHI());
      FLOAT P_X2 = PMU2.PT()*TMATH::COS(PMU2.PHI());
      FLOAT P_X3 = PMU3.PT()*TMATH::COS(PMU3.PHI());

      FLOAT P_Y1 = PMU1.PT()*TMATH::SIN(PMU1.PHI());
      FLOAT P_Y2 = PMU2.PT()*TMATH::SIN(PMU2.PHI());
      FLOAT P_Y3 = PMU3.PT()*TMATH::SIN(PMU3.PHI());

      FLOAT ET_X = MET.PT()*TMATH::COS(MET.PHI());
      FLOAT ET_Y = MET.PT()*TMATH::SIN(MET.PHI());

      FLOAT MTWZ = SQRT(POW((PMU1.PT() + PMU2.PT() + PMU3.PT() + MET.E()),2)  - POW((P_X1 + P_X2 + P_X3 + ET_X ),2) -  POW((P_Y1 + P_Y2 + P_Y3 + ET_Y ),2));

      IF ((MET_ET/1000 >30.&& MET_ET/1000. < 80.) && (mtwz > 120 && mtwz <350)){
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
          tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363491.lllv.2lep.root"); // WZ  MC_actual
      //   tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_363490.llll.2lep.root");  // ZZ MC
      //  tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_3641*.Zmumu_PTV0_70_CVetoBVeto.2lep.root");  //  Z + jet MC  
      //    tchain->Add("/cephfs/user/etoerne/ATLASOpen13TevData/2lep/MC/mc_392217.C1N2_WZ_400p0_0p0_3L_2L7.2lep.root"); // M1 
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
