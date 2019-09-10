
//  ---------------------------------------------------------------------------------
//
//   -- General macro for the analysis of radioactive decays using silicon detectors
//      Author: Jose Luis Rodriguez Sanchez, Daria Kostyleva
//      Last Update: 12/07/19
//      Comments: tracking in-flight decay products (30Cl->29S+p for instanse) using 
//      silicon microstrip AMS detectors
	
//  ---------------------------------------------------------------------------------
//   Usage: 
//      > root -l ana_si_detectors.C
//  ---------------------------------------------------------------------------------


#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
using namespace std;


void ana_si_detectors_kinem(){

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();
    // ------------------------------------------------------------------------

    //DEBUG  (optional)   -----------------------------------------------------
    gDebug = 0;

    //STYLE   -----------------------------------------------------------------		
    //gROOT->SetStyle("Default");
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);

    //DEFINE THE INPUT FILE  --------------------------------------------------

    //TString file_in = "/data.local2/G4_sim_momenta/sim_out_500k_110719.root";
    //TString file_in = "./sim_out.root";
    TString file_in = "./sim_out_kinem.root";
    TFile *file0 = TFile::Open(file_in);
    TTree* Tree0 = (TTree*)file0->Get("evt");
    Long64_t nevents = Tree0->GetEntries();
    std::cout<<"Number of entries: "<<nevents<<std::endl;

    TString file_hist_out = "./root_with_hist_kinem/out_hist.root";
    TFile * out_hist = new TFile(file_hist_out,"RECREATE");
    //TFile * out_hist = new TFile("out_hist_1m.root","RECREATE");

    //HISTOGRAMS DEFINITION

    TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.76,36.05);
   //   TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.6,36.2);
    TH1F* h_p_pz = new TH1F("h_p_pz","Longitudinal momentum of proton (GeV/c)",500,1.1,1.4);
    TH1F* h_s_pt = new TH1F("h_s_pt","Transverse momentum of HI (GeV/c)",500,0.0,0.08);
    TH1F* h_p_pt = new TH1F("h_p_pt","Transverse momentum of proton (GeV/c)",500,0.0,0.08);
    TH1F* h_s_px = new TH1F("h_s_px","Px of HI (GeV/c)",500,-0.1,0.1);
    TH1F* h_p_px = new TH1F("h_p_px","Px of proton (GeV/c)",500,-0.1,0.1);
    TH1F* ang_s_p = new TH1F("ang_s_p","Angle between 29S and proton (mrad)",1000,47,51);
    TH1F* ang_s_p_cut = new TH1F("ang_s_p_cut","Angle between 29S and proton (mrad) with selecton on pz s",500,40,55);
    TH2F* ang_momz = new TH2F("ang_momz","Correlation between long mom of HI and angle HI-p",500,35.76,36.05,500,0,55);
    TH2F* ang_momt = new TH2F("ang_momt","Correlation between transv mom of HI and angle HI-p",500,0.0,0.08,500,0,55);
    TH2F* ang_momtp = new TH2F("ang_momtp","Correlation between transv mom of proton and angle HI-p",500,0.0,0.08,500,0,55);
    TH2F* h_pt_pt = new TH2F("h_pt_pt","Correlation between transverse momenta of products",500,0.0,0.08,500,0.0,0.08);
    TH2F* h_pz_pz = new TH2F("h_pz_pz","Correlation between longitudinal momenta of products",500,1.1,1.4,500,35.78,36);
    TH2F* h_pt_pz = new TH2F("h_pt_pz","Correlation between long mom of HI and transv mom of proton",500,0.0,0.08,500,35.78,36);
    TH2F* h_pz_pt = new TH2F("h_pz_pt","Correlation between transv mom of HI and long mom of proton",500,1.1,1.4,500,0.0,0.08);
    TAxis* xaxis = new TAxis();
    TH2F* h_corr = new TH2F("h_corr","Correlation plot reflecting entries",500,35.76,36.05,500,1,500);
    TH2F* h_pz_pt_s = new TH2F("h_pz_pt_s","pz 29S - x, pt 29S - y",500,35.76,36.05,500,0.0,0.08);
    //energies
    TH1F* h_s_en = new TH1F("h_s_en","Total energy of 29S (GeV)",500,44.5,45.5);
    TH1F* h_p_en = new TH1F("h_p_en","Total energy of proton (GeV)",500,1.,2.0);
    TH1F* h_s_kin = new TH1F("h_s_kin","Kinetic energy of 29S (GeV)",500,17.5,18.5);
    TH1F* h_p_kin = new TH1F("h_p_kin","Kinetic energy of proton (GeV)",500,0.,1.0);
    TH1F* h_s_de = new TH1F("h_s_de","Energy loss of 29S (GeV)",500,0.0,0.001);
    TH1F* h_p_de = new TH1F("h_p_de","Energy loss of proton (GeV)",500,0.0,0.0001);

    TH1F* h_qval = new TH1F("h_qval","Q-value via kin energy difference (MeV)",500,1.4,2.4);
    TH1F* h_qval1 = new TH1F("h_qval1","Q-value via kin energy difference (MeV)",500,1.4,2.4);
    TH1F* h_qval2 = new TH1F("h_qval2","Q-value via proton pcm (MeV)",500,1.98,2.02);
    TH1F* h_qval4 = new TH1F("h_qval4","Q-value via proton pcm (MeV)",500,1.98,2.02);
    TH1F* h_qval3 = new TH1F("h_qval3","Q-value via HI pcm (MeV)",500,1.8,2.2);
    TH1F* h_qval5 = new TH1F("h_qval5","Q-value via HI pcm (MeV)",500,1.8,2.2);
    TH1F* h_qval6 = new TH1F("h_qval6","Q-value via HI pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval7 = new TH1F("h_qval7","Q-value via HI pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval8 = new TH1F("h_qval8","Q-value via proton pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval9 = new TH1F("h_qval9","Q-value via proton pcm at FRS (MeV)",500,-0.1,3.);

    //Monte-Carlo Track (input)
    TClonesArray* MCTrackCA;  
    R3BMCTrack** track;
    MCTrackCA = new TClonesArray("R3BMCTrack",5);
    TBranch *branchMCTrack = Tree0->GetBranch("MCTrack");
    branchMCTrack->SetAddress(&MCTrackCA);

    //Tracker Hits (input)
    TClonesArray* TraCA;   
    R3BTraPoint** Tra;
    TraCA = new TClonesArray("R3BTraPoint",5);
    TBranch *branchTra = Tree0->GetBranch("TraPoint");
    branchTra->SetAddress(&TraCA);

    Int_t primary=0;
    TVector3 momentum;
    TVector3 momentum_s;
    TVector3 momentum_p;
    Int_t traPerEvent=0;
    Int_t MCtracksPerEvent=0;

    //***** Variables for tracking ********
    Double_t x[2][4], y[2][4], z[2][4]; //1st index - particle (0 - 29S, 1 - proton), 2nd - detector number (there are 4 silicons)
    Double_t dx[2], dy[2], dz[2], rad[2];
    Double_t cos_ang, ang_p1S, ang_p1S_max = 0.0;
    Double_t s_pz, p_pz, p_pt, s_pt, s_px, s_py, p_px, p_py;
    Double_t det01, dist01, de_p = 0., de_s = 0.;
    Double_t s_pz_max = 35.9, s_pz_min = 35.9;
    

    //***** Values for cheking qvalue reproducement, taken from asciigenerator *******//
    Double_t q_kin = 0.;   // I calculate by the formula using difference of kinetic energies
    Double_t qvalue_in = 0.002; // I set initially
    Double_t Mp = 0.938272297;
    Double_t Mhi = 29.*0.93149406-0.003156;
    Double_t Pmom, Phi, Pp;
    Double_t Mmom = Mhi + Mp + qvalue_in;


    //***** Calculate qvalue from obtained plab of HI (or proton) and converted to pcm and put into pcm formula
    Double_t q_pcm_p=0. , q_pcm_s = 0.; 
    Double_t g = 1.66332, b = -0.799092, pcm_s = 0., pcm_p = 0.; // beta and gamma from simulation directly
    Double_t e_hi_lab = 0.; 
    Double_t e_p_lab = 0.; 
    Double_t s_pz_1, p_pz_1; //momenta components of hi and p for reverce lorenty boost
    Double_t s_pz_1_ecorr, p_pz_1_ecorr, q_pcm_s_ecorr, q_pcm_p_ecorr, pcm_s_ecorr, pcm_p_ecorr;
    Double_t pcm_s_frs = 0., q_pcm_s_frs = 0.,pcm_s_frs1 = 0., q_pcm_s_frs1 = 0.  ;
    Double_t pcm_p_frs = 0., q_pcm_p_frs = 0.,pcm_p_frs1 = 0., q_pcm_p_frs1 = 0.  ;
      //******* Defining angles for isotropic distribution ************
    //Double_t costheta = 2.*gRandom->Uniform(0,1)-1.0;
    //Double_t sintheta = std::sqrt((1.0 - costheta)*(1.0 + costheta));
    //Double_t phi  = 2*TMath::Pi()*gRandom->Uniform(0,1);

    //****** Conditions on the event selection *****
    Double_t min_angle = 2.0;   // minimum angle between trajectories is 2 mrad
    Double_t same_strip = 0.03; // minimum distance 300mum between coordinates of hi and p, if less than they hit the same strip 
    Double_t max_dist = 0.018;  // maximum distance between trajectories 180 mum, if greater hi and p don't come from the same mother

    //************* MAIN LOOP OVER EVENTS *************************************
    //for(Int_t i=0;i<100;i++){
    for(Int_t i=0;i<nevents;i++){
       cout << "Event " << i << endl;
       if(i%10000 == 0) printf("Event:%i\n",i);

       MCTrackCA->Clear();
       TraCA->Clear();

       Tree0->GetEvent(i);
       MCtracksPerEvent = MCTrackCA->GetEntries();
       traPerEvent = TraCA->GetEntries();

       if(traPerEvent>0) {
	        Tra = new R3BTraPoint*[traPerEvent];
	        for(Int_t j=0;j<traPerEvent;j++){
	           Tra[j] = new R3BTraPoint;
             Tra[j] = (R3BTraPoint*) TraCA->At(j);      
	        }
       }
		
       if(MCtracksPerEvent>0) {
	        track = new R3BMCTrack*[MCtracksPerEvent];
	        for(Int_t j=0;j<MCtracksPerEvent;j++){
       	     track[j] = new R3BMCTrack;
	           track[j] = (R3BMCTrack*) MCTrackCA->At(j);      
	        }
       }

       // loop in MC mother tracks per event. Now we have 29S and proton here
       // because I generate decay myself, I have to register daughters as primary ions, because there is no p-emission in G4
       // but the masses which are attached to them are wrong! just p_mass*A. I don't use those masses, I have to re-assign the real masses each time

       for(Int_t h=0;h<MCtracksPerEvent;h++){
         // cout << "MCTrack " <<h<<" in event "<<i<<endl; 
		      if(track[h]->GetMotherId()<0){
			       track[h]->GetMomentum(momentum);	
		         //h1_T->Fill(track[h]->GetPdgCode()-1000000000);
             //h2_T->Fill(track[h]->GetEnergy()*1000.-track[h]->GetMass()*1000.);
			       primary++;	
             if(track[h]->GetPdgCode()==1000160290) {
              // !!!!!!!!!!!!!!!!!!GetMass here for HI is wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             }
             if(track[h]->GetPdgCode()==2212) {
              //this is proton, so assigned mass by G4 is correct, but I don't use it anyway
             }            
          }   
        }//for MC

       // *************** LOOP over tracks per event! ****************************************
       // *************** GetDetCopyID() - detector number, 1,2,3 or 4
       // here we always (almost) have 8 tracks per events because our 2 particles pass trough all 4 detectors 2x4=8
       // in kinematic case I know that I generate 2 particles and thus I don't look for coincidences
       for(Int_t h=0;h<traPerEvent;h++){          
           //******* 29S *******************
           if(Tra[h]->GetPdi()==1000160290 && Tra[h]->GetDetCopyID()==1){
              Tra[h]->MomentumOut(momentum_s);
              s_pt = sqrt(Tra[h]->GetPxOut()*Tra[h]->GetPxOut() + Tra[h]->GetPyOut()*Tra[h]->GetPyOut());
              s_px = Tra[h]->GetPxOut();
              s_py = Tra[h]->GetPyOut();
              s_pz = Tra[h]->GetPzOut();
              x[0][0] = Tra[h]->GetXOut();
              y[0][0] = Tra[h]->GetYOut();
              z[0][0] = Tra[h]->GetZOut();
              de_s = Tra[h]->GetEloss();
              e_hi_lab = sqrt(s_px*s_px + s_py*s_py + s_pz*s_pz + Mhi*Mhi);
             // printf("e_hi_lab = %f\n",e_hi_lab);
            //  cout << "Detector " << Tra[h]->GetDetCopyID() << endl;
              if(s_pz > s_pz_max) s_pz_max = s_pz;
              if(s_pz < s_pz_min && s_pz > 0.) s_pz_min = s_pz;
             // printf("pz is %f, pz max is %f, pz min is %f\n",s_pz, s_pz_max,s_pz_min);
           }
           //******* proton ****************
           if(Tra[h]->GetPdi()==2212 && Tra[h]->GetDetCopyID()==1){
              Tra[h]->MomentumOut(momentum_p);
              p_pt = sqrt(Tra[h]->GetPxOut()*Tra[h]->GetPxOut() + Tra[h]->GetPyOut()*Tra[h]->GetPyOut());
              p_pz = Tra[h]->GetPzOut();
              p_px = Tra[h]->GetPxOut();
              p_py = Tra[h]->GetPyOut();            
              x[1][0] = Tra[h]->GetXOut();
              y[1][0] = Tra[h]->GetYOut();
              z[1][0] = Tra[h]->GetZOut(); 
              de_p = Tra[h]->GetEloss(); 
              e_p_lab = sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz + Mp*Mp);     
              //printf("e_p_lab = %f\n",e_p_lab);       
            }
       }//for tracker

       cout << "Tracks per event " << traPerEvent << endl;
       // ************** END OF LOOP over tracks per event ************************************

       //*********** Calculation of the angles between tracks *************
       //*********** For only 1st silicon detector now, i.e. second index 0 ****

       rad[0] = sqrt(x[0][0]*x[0][0] + y[0][0]*y[0][0] + z[0][0]*z[0][0]);
       rad[1] = sqrt(x[1][0]*x[1][0] + y[1][0]*y[1][0] + z[1][0]*z[1][0]);

       // cos of angles = scalar product/their lengths
       cos_ang = (x[0][0]*x[1][0] + y[0][0]*y[1][0] + z[0][0]*z[1][0])/(rad[0]*rad[1]);
       ang_p1S = acos(cos_ang)*1000.;
       if(ang_p1S > ang_p1S_max) ang_p1S_max = ang_p1S;

       if(ang_p1S < min_angle) continue; //angle less than 2 mrad
       //if decay products fall into the same strip i.e. dead zones
       if(abs(x[0][0]-x[1][0]) < same_strip) continue; //300 mum 
       if(abs(y[0][0]-y[1][0]) < same_strip) continue;

       //minimum distance between tracks, 0 - 29S, 1 - proton should not exceed 180 um
       det01 = (x[0][0]-x[1][0])*(y[0][0]*z[1][0]-y[1][0]*z[0][0]) + (y[0][0]-y[1][0])*(-x[0][0]*z[1][0]+x[1][0]*z[0][0]) + (z[0][0]-z[1][0])*(x[0][0]*y[1][0]-x[1][0]*z[0][0]);
       //det12=(x(1,1)-x(1,2))*(dy(1)*dz(2)-dy(2)*dz(1)) +(y(1,1)-y(1,2))*(-dx(1)*dz(2)+dx(2)*dz(1)) +(z(1,1)-z(1,2))*(dx(1)*dy(2)-dx(2)*dz(1))
       dist01 = abs(det01/(rad[0]*rad[1]*sqrt(1.0-cos_ang*cos_ang)));
       if(dist01 > max_dist) continue;
       //******************************************************************

       //***** cheking q-value reproducement ***********************************************
       // **** q-value via kinetic energy difference ***************************************
       Pmom = sqrt(0.618*30*(0.618*30+ 2*(Mhi + Mp + qvalue_in)));
       Phi = sqrt(s_px*s_px + s_py*s_py + s_pz*s_pz); //3-momentum in lab
       Pp = sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz); //3-momentum in lab
       q_kin = sqrt(Pp*Pp+Mp*Mp) + sqrt(Phi*Phi + Mhi*Mhi) - sqrt(Pmom*Pmom + Mmom*Mmom) - (Mp + Mhi -Mmom);
       //cout << "Qvalue = " <<q_kin << endl;
       //***********************************************************************************


       // **** q-value from pcm of HI and proton separately ********************************
       //calculate q-value backwards, i.e. from plab of hi.

       //first Lorentz-boost to cm, momenta are in cartesian coordinates
       s_pz_1 = (s_pz + g*b*((g*b*s_pz)/(g+1) + e_hi_lab));
       p_pz_1 = (p_pz + g*b*((g*b*p_pz)/(g+1) + e_p_lab));
       //cout << "s_pz_1 = " << s_pz_1 << endl;
       //cout << "p_pz_1 = " << p_pz_1 << endl;

       //from cartezian coordinates to spherical
       pcm_s = sqrt(s_px*s_px + s_py*s_py + s_pz_1*s_pz_1); //3-momentum in cm
       pcm_p = sqrt(p_px*p_px + p_py*p_py + p_pz_1*p_pz_1); //3-momentum in cm, these two are not exactly equal because no energy loss is taken into account
       //cout << "pcm_s = " << pcm_s << endl;
       //cout << "pcm_p = " << pcm_p << endl;
       q_pcm_s = 4*Mmom*Mmom*pcm_s*pcm_s/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       q_pcm_p = 4*Mmom*Mmom*pcm_p*pcm_p/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));

       // energy-loss correction
       s_pz_1_ecorr = (s_pz + g*b*((g*b*s_pz)/(g+1) + e_hi_lab + de_s ));
       p_pz_1_ecorr = (p_pz + g*b*((g*b*p_pz)/(g+1) + e_p_lab  + de_p ));
       pcm_s_ecorr = sqrt(s_px*s_px + s_py*s_py + s_pz_1_ecorr*s_pz_1_ecorr); //3-momentum in cm
       pcm_p_ecorr = sqrt(p_px*p_px + p_py*p_py + p_pz_1_ecorr*p_pz_1_ecorr);
       q_pcm_s_ecorr = 4*Mmom*Mmom*pcm_s_ecorr*pcm_s_ecorr/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       q_pcm_p_ecorr = 4*Mmom*Mmom*pcm_p_ecorr*pcm_p_ecorr/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       
       // in experiment at FRS only long mom of HI is measured
       //cout << "s_px " << s_px << endl;
       //cout << "s_py " << s_py << endl;
       //cout << "s_pz_1 " << s_pz_1 << endl;
       //cout << "p_pz_1 " << p_pz_1 << endl;
       //pcm_s_frs = s_pz_1/costheta;
       //****** HI ***************
       pcm_s_frs = sqrt(s_pz_1*s_pz_1);
       pcm_s_frs1 = sqrt(s_px*s_px + s_py*s_py);
       //****** Proton *************
       pcm_p_frs = sqrt(p_pz_1*p_pz_1);
       pcm_p_frs1 = sqrt(p_px*p_px + p_py*p_py);
       //cout << "pcm_s_frs  " << pcm_s_frs  << endl;

       //****** HI ***************
       q_pcm_s_frs = 4*Mmom*Mmom*pcm_s_frs*pcm_s_frs/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       q_pcm_s_frs1 = 4*Mmom*Mmom*pcm_s_frs1*pcm_s_frs1/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));

       //****** Proton *************
       q_pcm_p_frs = 4*Mmom*Mmom*pcm_p_frs*pcm_p_frs/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       q_pcm_p_frs1 = 4*Mmom*Mmom*pcm_p_frs1*pcm_p_frs1/((Mmom + Mp - Mhi)*(Mmom - Mp + Mhi)*(Mmom + Mp + Mhi));
       //**********************************************************************************



       //*********** Filling histos ***************************************
       if(s_pz > 35.8755 && s_pz < 35.94) ang_s_p_cut->Fill(ang_p1S);
       ang_s_p->Fill(ang_p1S);
       //if(s_pz > 35.8755 && s_pz < 35.94) ang_momz->Fill(s_pz,ang_p1S);
       ang_momz->Fill(s_pz,ang_p1S);
       ang_momt->Fill(s_pt,ang_p1S);
       ang_momtp->Fill(p_pt,ang_p1S);
       h_s_pz->Fill(s_pz);
       h_p_pz->Fill(p_pz);
       h_s_pt->Fill(s_pt);
       h_p_pt->Fill(p_pt);
       h_s_px->Fill(s_px);
       h_p_px->Fill(p_px);
       h_pt_pt->Fill(p_pt,s_pt);
       h_pz_pz->Fill(p_pz,s_pz);
       h_pt_pz->Fill(p_pt,s_pz);
       h_pz_pt->Fill(p_pz,s_pt);
       h_pz_pt_s->Fill(s_pz,s_pt);
       h_s_de->Fill(de_s);
       h_p_de->Fill(de_p); 
       h_qval->Fill(q_kin*1000);
       h_qval1->Fill(q_kin*1000 + de_p*1000 + de_s*1000);
       h_qval2->Fill(q_pcm_p*1000);
       h_qval3->Fill(q_pcm_s*1000);
       h_qval4->Fill(q_pcm_p_ecorr*1000);
       h_qval5->Fill(q_pcm_s_ecorr*1000);
       h_qval6->Fill(q_pcm_s_frs*1000);
       h_qval7->Fill(q_pcm_s_frs1*1000);
       h_qval8->Fill(q_pcm_p_frs*1000);
       h_qval9->Fill(q_pcm_p_frs1*1000);
       h_s_en->Fill(e_hi_lab);
       h_p_en->Fill(e_p_lab);
       h_s_kin->Fill(e_hi_lab - Mhi);
       h_p_kin->Fill(e_p_lab - Mp);


       //******************************************************************
       //cout << "Num of prim "<<primary << endl;

       if(traPerEvent>0)         delete[] Tra;
       if(MCtracksPerEvent>0)    delete[] track;

    }//loop for nevents


    xaxis = h_s_pz->GetXaxis();
    for(Int_t i=0; i<500; i++){
      for(Int_t j=0; j<500; j++){
        h_corr->Fill(xaxis->GetBinCenter(i),ang_momz->GetBinContent(i,j));
      }
    }
    printf("max angle betweeen HI and protons is %f mrad \n", ang_p1S_max);
    printf("pz max is %f GeV/c, pz min is %f GeV/c \n",s_pz_max,s_pz_min);
    printf("delta p of HI = %f GeV/c \n",s_pz_max-s_pz_min);
    printf("FRS acceptance = %f\n",h_s_pz->GetMean()/100*1); // FRS acceptance is 2% or +-1%
    printf("Super-FRS acceptance = %f\n",h_s_pz->GetMean()/100*2.5); // +-2.5%


    //MC TRACK CANVAS
    TCanvas* c1 = new TCanvas("MCTrack","MCTrack",0,0,720,900);
    c1->SetFillColor(0);
    c1->SetFrameFillColor(0);
    c1->cd();
    c1->Divide(1,2);
    c1->cd(1);
    ang_momtp->Draw();
    c1->cd(2);
    h_pz_pt_s->Draw();

    TCanvas* c2 = new TCanvas("angles","angles",0,0,720,900);
    ang_s_p->Draw();
    ang_s_p_cut->SetLineColor(2);
    ang_s_p_cut->Draw("same");

  
  //in-flight decay tracking pics
    TCanvas* c3 = new TCanvas("momenta1", "Momenta in lab",0,0,800,900);
    c3->Divide(2,3);
    c3->cd(1);  
    h_s_pz->GetXaxis()->SetTitle("Pz of HI (GeV/c)");
  //  h_s_pz->GetYaxis()->SetTitle("Probability from KS");
    h_s_pz->SetFillColor(kBlue-8);
    h_s_pz->Draw("bar");
    c3->cd(2);
    h_p_pz->GetXaxis()->SetTitle("Pz of proton (GeV/c)");
    h_p_pz->SetFillColor(kBlue-8);
    h_p_pz->Draw("bar");
    c3->cd(3);
    h_s_px->GetXaxis()->SetTitle("Px of HI (GeV/c)");
    h_s_px->Draw();
   
    c3->cd(4);
    h_p_px->GetXaxis()->SetTitle("Px of proton (GeV/c)");
    h_p_px->Draw();
 

    c3->cd(5);
    h_s_pt->GetXaxis()->SetTitle("Pt of HI (GeV/c)");
    h_s_pt->Draw();
    c3->cd(6);
    h_p_pt->GetXaxis()->SetTitle("Pt of proton (GeV/c)");
    h_p_pt->Draw();

    c3->cd(1);
    TH1D *project1 = ang_momz->ProjectionX();
    project1->SetLineColor(2);
    //project1->Draw("same");

    TCanvas * c4 = new TCanvas("momenta2", "Momenta in lab",800,0,800,900);
    c4->Divide(2,2);
    c4->cd(1);
    h_pt_pt->GetXaxis()->SetTitle("Pt of proton (GeV/c)");
    h_pt_pt->GetYaxis()->SetTitle("Pt of HI (GeV/c)");
    h_pt_pt->Draw("colz");
    c4->cd(2);
    h_pz_pz->GetXaxis()->SetTitle("Pz of proton (GeV/c)");
    h_pz_pz->GetYaxis()->SetTitle("Pz of HI (GeV/c)");
    h_pz_pz->Draw("colz");
    c4->cd(3);
    h_pt_pz->GetXaxis()->SetTitle("Pz of proton (GeV/c)");
    h_pt_pz->GetYaxis()->SetTitle("Pt of HI (GeV/c)");
    h_pt_pz->Draw("colz");
    c4->cd(4);
    h_pz_pt->GetXaxis()->SetTitle("Pt of proton (GeV/c)");
    h_pz_pt->GetYaxis()->SetTitle("Pz of HI (GeV/c)");
    h_pz_pt->Draw("colz");

    TCanvas * c5 = new TCanvas("momenta3", "Momenta in lab",1600,0,800,900);
    c5->Divide(2,2);
    c5->cd(1);
    ang_s_p->GetXaxis()->SetTitle("Angle between HI and proton in lab (mrad)");
    ang_s_p->Draw();
    //ang_s_p_cut->SetLineColor(2);
    //ang_s_p_cut->Draw("same");
    c5->cd(2);
    ang_momz->GetXaxis()->SetTitle("Pz of HI (GeV/c)");
    ang_momz->GetYaxis()->SetTitle("Angle between HI and proton in lab (mrad)");
    ang_momz->Draw("colz");
    c5->cd(3);
    ang_momt->GetXaxis()->SetTitle("Pt of HI (GeV/c)");
    ang_momt->GetYaxis()->SetTitle("Angle between HI and proton in lab (mrad)");
    ang_momt->Draw("colz");
    c5->cd(4);
    h_corr->GetXaxis()->SetTitle("Pz of HI (GeV/c)");
    h_corr->GetYaxis()->SetTitle("Bin content of correlation plot: pz pf HI vs angle");
    h_corr->Draw("col");
 

    TCanvas * c6 = new TCanvas("energy", "Enenrgy in lab",0,900,800,900);
    c6->Divide(2,3);
    c6->cd(1);
    h_s_en->GetXaxis()->SetTitle("Total energy of HI (GeV)");
    h_s_en->Draw();
    c6->cd(2);
    h_p_en->GetXaxis()->SetTitle("Total energy of proton (GeV)");
    h_p_en->Draw();   
    c6->cd(3);
    h_s_kin->GetXaxis()->SetTitle("Kinetic energy of HI (GeV)");
    h_s_kin->Draw();
    c6->cd(4);
    h_p_kin->GetXaxis()->SetTitle("Kinetic energy of proton (GeV)");
    h_p_kin->Draw();
    c6->cd(5);
    h_p_de->GetXaxis()->SetTitle("Energy loss of proton (GeV)");
    h_p_de->Draw();
    c6->cd(6);
   // h_qval->Draw();
   // h_qval1->SetLineColor(2);
  //  h_qval1->Draw("same");
    h_s_de->GetXaxis()->SetTitle("Energy loss of HI (GeV)");
    h_s_de->Draw();


    TCanvas * c7 = new TCanvas("qvalues", "Decay energies calculated differently",900,900,800,900);
    gStyle->SetOptFit();
    c7->Divide(1,5);
    c7->cd(1);
    h_qval->GetXaxis()->SetTitle("Qvalue calculated via kinetic energies after tracking (MeV)");
    h_qval->GetXaxis()->SetLabelSize(0.05);
    h_qval->GetXaxis()->SetTitleSize(0.05);
    h_qval1->GetXaxis()->SetTitle("Qvalue calculated via kinetic energies after tracking (MeV)");
    h_qval1->GetXaxis()->SetLabelSize(0.05);
    h_qval1->GetXaxis()->SetTitleSize(0.05);
    h_qval1->SetLineColor(8);
    h_qval1->Draw();
    h_qval->Draw("same");
    TLegend *l7_1 = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l7_1->AddEntry("h_qval","Q-value before energy loss correction","l");
    l7_1->AddEntry("h_qval1","Q-value after energy loss correction","l");
    l7_1->Draw("same");

    c7->cd(2);
    h_qval2->GetXaxis()->SetTitle("Qvalue calculated via pcm of proton (MeV)");
    h_qval2->GetXaxis()->SetLabelSize(0.05);
    h_qval2->GetXaxis()->SetTitleSize(0.05);
    h_qval2->Draw();
    h_qval4->SetLineColor(8);
    h_qval2->Fit("gaus","","",1.999,2.001);
    h_qval4->Draw("same");
    TLegend *l8_1 = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l8_1->AddEntry("h_qval2","Q-value before energy loss correction","l");
    l8_1->AddEntry("h_qval4","Q-value after energy loss correction","l");
    //l8_1->AddEntry("gaus", "Gaussian fit", "l");
    l8_1->Draw("same");

    c7->cd(3);
    h_qval3->GetXaxis()->SetTitle("Qvalue calculated via pcm of HI (MeV)");
    h_qval3->GetXaxis()->SetLabelSize(0.05);
    h_qval3->GetXaxis()->SetTitleSize(0.05);
    h_qval3->Draw();
    h_qval5->SetLineColor(8);
    h_qval3->Fit("gaus");
    h_qval5->Draw("same");
    TLegend *l8_2 = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l8_2->AddEntry("h_qval3","Q-value before energy loss correction","l");
    l8_2->AddEntry("h_qval5","Q-value after energy loss correction","l");
    //l8_2->AddEntry("gaus", "Gaussian fit", "l");
    l8_2->Draw("same");

    c7->cd(4);
    h_qval6->GetXaxis()->SetTitle("Qvalue calculated via pcm of HI with diff plab components (MeV)");
    h_qval6->GetXaxis()->SetLabelSize(0.05);
    h_qval6->GetXaxis()->SetTitleSize(0.05);
    h_qval6->Draw();
    h_qval6->SetLineColor(6);
    h_qval7->Draw("sames");
    h_qval7->Fit("gaus","","", 1.95, 2.1);
    TLegend *l7 = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l7->AddEntry("h_qval6","Q-value via p longitudinal of HI","l");
    l7->AddEntry("h_qval7","Q-value via p transverse of HI","l");
    //l7->AddEntry("gaus", "Gaussian fit", "l");
    l7->Draw("same");
    //c7->SaveAs("canvas_with_qvalue.root");

    c7->cd(5);
    h_qval8->GetXaxis()->SetTitle("Qvalue calculated via pcm of proton with diff plab components (MeV)");
    h_qval8->GetXaxis()->SetLabelSize(0.05);
    h_qval8->GetXaxis()->SetTitleSize(0.05);
    h_qval8->Draw();
    h_qval8->SetLineColor(6);
    h_qval9->Draw("sames");
   // h_qval9->Fit("gaus","","", 1.95, 2.1);
   // TLegend *l7 = new TLegend(0.1,0.7,0.48,0.9);


    // writing histograms into root file 
    ang_s_p->Write();
    ang_momz->Write();
    ang_momt->Write();
    ang_momtp->Write();
    h_s_pz->Write();
    h_p_pz->Write();
    h_s_pt->Write();
    h_p_pt->Write();
    h_s_px->Write();
    h_p_px->Write();
    h_pt_pt->Write();
    h_pz_pz->Write();
    h_pt_pz->Write();
    h_pz_pt->Write();
    h_pz_pt_s->Write();
    h_corr->Write();
    h_p_de->Write();
    h_s_de->Write();
    h_qval->Write();
    h_qval1->Write();
    h_qval2->Write();
    h_qval3->Write();
    h_qval4->Write();
    h_qval5->Write();
    h_qval6->Write();
    h_qval7->Write();
    h_s_en->Write();
    h_p_en->Write();
    h_s_kin->Write();
    h_p_kin->Write();
    c1->Write();c2->Write();c3->Write();c4->Write();c5->Write();c6->Write();c7->Write();



    // -----   Finish   -------------------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    cout << endl << endl;
    cout << "Macro finished succesfully." << endl;
    cout << "Real time " << rtime << " s, CPU time " << ctime << " s" << endl;
    cout << endl;
    // ------------------------------------------------------------------------
}

