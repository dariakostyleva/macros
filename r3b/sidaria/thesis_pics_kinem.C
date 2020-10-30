
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
#include "TLatex.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
using namespace std;


void thesis_pics_kinem(){

    // -----   Timer   --------------------------------------------------------
    TStopwatch timer;
    timer.Start();
    // ------------------------------------------------------------------------

    //DEBUG  (optional)   -----------------------------------------------------
    gDebug = 0;

    //STYLE   -----------------------------------------------------------------		
    //gROOT->SetStyle("Default");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    //DEFINE THE INPUT FILE  --------------------------------------------------

    TString file_in = "/data.local2/G4_sim_momenta/kinem_case/sim_out_500k_110719.root";
    //TString file_in = "./sim_out.root";
    //TString file_in = "./sim_out_kinem.root";
    TFile *file0 = TFile::Open(file_in);
    TTree* Tree0 = (TTree*)file0->Get("evt");
    Long64_t nevents = Tree0->GetEntries();
    std::cout<<"Number of entries: "<<nevents<<std::endl;

   // TString file_hist_out = "./root_with_hist_kinem/out_hist_500k_110520_thesis.root";
   // TFile * out_hist = new TFile(file_hist_out,"RECREATE");
    //TFile * out_hist = new TFile("out_hist_1m.root","RECREATE");

    //HISTOGRAMS DEFINITION

    TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.76,36.05);
   //   TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.6,36.2);
   // TH1F* h_p_pz = new TH1F("h_p_pz","Longitudinal momentum of proton (GeV/c)",500,1.1,1.4);
   // TH1F* h_s_pt = new TH1F("h_s_pt","Transverse momentum of HI (GeV/c)",500,0.0,0.08);
   // TH1F* h_p_pt = new TH1F("h_p_pt","Transverse momentum of proton (GeV/c)",500,0.0,0.08);
   // TH1F* h_s_px = new TH1F("h_s_px","Px of HI (GeV/c)",500,-0.1,0.1);
   // TH1F* h_p_px = new TH1F("h_p_px","Px of proton (GeV/c)",500,-0.1,0.1);
    TH1F* ang_s_p = new TH1F("ang_s_p","Angle between 29S and proton (mrad)",300,0,60);
    //TH1F* ang_s_p_cut = new TH1F("ang_s_p_cut","Angle between 29S and proton (mrad) with selecton on pz s",500,40,55);
    TH2F* ang_momz = new TH2F("ang_momz","Correlation between long mom of HI and angle HI-p",300,0,55,600,35.76,36.05);
    //TH2F* ang_momt = new TH2F("ang_momt","Correlation between transv mom of HI and angle HI-p",500,0.0,0.08,500,0,55);
    //TH2F* ang_momtp = new TH2F("ang_momtp","Correlation between transv mom of proton and angle HI-p",500,0.0,0.08,500,0,55);
   // TH2F* h_pt_pt = new TH2F("h_pt_pt","Correlation between transverse momenta of products",500,0.0,0.08,500,0.0,0.08);
   // TH2F* h_pz_pz = new TH2F("h_pz_pz","Correlation between longitudinal momenta of products",500,1.1,1.4,500,35.78,36);
   // TH2F* h_pt_pz = new TH2F("h_pt_pz","Correlation between long mom of HI and transv mom of proton",500,0.0,0.08,500,35.78,36);
   // TH2F* h_pz_pt = new TH2F("h_pz_pt","Correlation between transv mom of HI and long mom of proton",500,1.1,1.4,500,0.0,0.08);
   // TAxis* xaxis = new TAxis();
   // TH2F* h_corr = new TH2F("h_corr","Correlation plot reflecting entries",500,35.76,36.05,500,1,500);
   // TH2F* h_pz_pt_s = new TH2F("h_pz_pt_s","pz 29S - x, pt 29S - y",500,35.76,36.05,500,0.0,0.08);
    //energies
   /* TH1F* h_s_en = new TH1F("h_s_en","Total energy of 29S (GeV)",500,44.5,45.5);
    TH1F* h_p_en = new TH1F("h_p_en","Total energy of proton (GeV)",500,1.,2.0);
    TH1F* h_s_kin = new TH1F("h_s_kin","Kinetic energy of 29S (GeV)",500,17.5,18.5);
    TH1F* h_p_kin = new TH1F("h_p_kin","Kinetic energy of proton (GeV)",500,0.,1.0);
    TH1F* h_s_de = new TH1F("h_s_de","Energy loss of 29S (GeV)",500,0.0,0.001);
    TH1F* h_p_de = new TH1F("h_p_de","Energy loss of proton (GeV)",500,0.0,0.0001);
    */
   // TH1F* h_qval = new TH1F("h_qval","Q-value via kin energy difference (MeV)",500,1.4,2.4);
   // TH1F* h_qval1 = new TH1F("h_qval1","Q-value via kin energy difference (MeV)",500,1.4,2.4);
    TH1F* h_qval2 = new TH1F("h_qval2","Q-value via proton pcm (MeV)",500,1.98,2.02);
   /* TH1F* h_qval4 = new TH1F("h_qval4","Q-value via proton pcm (MeV)",500,1.98,2.02);
    TH1F* h_qval3 = new TH1F("h_qval3","Q-value via HI pcm (MeV)",500,1.8,2.2);
    TH1F* h_qval5 = new TH1F("h_qval5","Q-value via HI pcm (MeV)",500,1.8,2.2);
    TH1F* h_qval6 = new TH1F("h_qval6","Q-value via HI pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval7 = new TH1F("h_qval7","Q-value via HI pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval8 = new TH1F("h_qval8","Q-value via proton pcm at FRS (MeV)",500,-0.1,3.);
    TH1F* h_qval9 = new TH1F("h_qval9","Q-value via proton pcm at FRS (MeV)",500,-0.1,3.);
    */
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

    bool hi_1, p_1, hi_2, p_2, hi_3, p_3, hi_4, p_4;  // flags to help looking for coincidences of hi and p in the 1st and 4th detectors 

    //************* MAIN LOOP OVER EVENTS *************************************
    //for(Int_t i=0;i<10000;i++){
    for(Int_t i=0;i<nevents;i++){
       //cout << "Event " << i << endl;
       hi_1 = false; p_1 = false;  // at the beggining of each event logic is reloaded
       hi_2 = false; p_2 = false;
       hi_3 = false; p_3 = false;
       hi_4 = false; p_4 = false;

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
              hi_1 = true;
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
            p_1 = true;
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
            //******* 29S in 4th silicon ****************
            if(Tra[h]->GetPdi()==1000160290 && Tra[h]->GetDetCopyID()==4) { 
              hi_4 = true; 
              x[0][3] = Tra[h]->GetXOut();
              y[0][3] = Tra[h]->GetYOut();
              z[0][3] = Tra[h]->GetZOut();
              //cout << " xyz for 29S in 4th detector " << x[0][3] << " " << y[0][3] << " " << z[0][3] <<endl;
            }

            //******* proton in 4th silicon ****************
            if(Tra[h]->GetPdi()==2212 && Tra[h]->GetDetCopyID()==4) { 
              p_4 = true; 
              x[1][3] = Tra[h]->GetXOut();
              y[1][3] = Tra[h]->GetYOut();
              z[1][3] = Tra[h]->GetZOut(); 
              //cout << " xyz for proton in 4th detector " << x[1][3] << " " << y[1][3] << " " << z[1][3] <<endl;
            }
       }//for tracker

       // ************** END OF LOOP over tracks per event ************************************

       //*********** Calculation of the angles between tracks *************
       //*********** between 1st and 4th silicon (second indexes 0 and 3)  ****


       // calculating coordinates of vectors belonging to hi (index 0) and p (index 1)
       dx[0] = x[0][3] - x[0][0]; dy[0] = y[0][3] - y[0][0]; dz[0] = z[0][3] - z[0][0];
       dx[1] = x[1][3] - x[1][0]; dy[1] = y[1][3] - y[1][0]; dz[1] = z[1][3] - z[1][0];

       rad[0] = sqrt(dx[0]*dx[0] + dy[0]*dy[0] + dz[0]*dz[0]);
       rad[1] = sqrt(dx[1]*dx[1] + dy[1]*dy[1] + dz[1]*dz[1]);
        //rad[0] = sqrt(x[0][0]*x[0][0] + y[0][0]*y[0][0] + z[0][0]*z[0][0]);
       //rad[1] = sqrt(x[1][0]*x[1][0] + y[1][0]*y[1][0] + z[1][0]*z[1][0]);

       // cos of angles = scalar product/their lengths
       cos_ang = (dx[0]*dx[1] + dy[0]*dy[1] + dz[0]*dz[1])/(rad[0]*rad[1]);
       //cos_ang = (x[0][0]*x[1][0] + y[0][0]*y[1][0] + z[0][0]*z[1][0])/(rad[0]*rad[1]);

       ang_p1S = acos(cos_ang)*1000.;
       if(ang_p1S > ang_p1S_max) ang_p1S_max = ang_p1S;

       //if(ang_p1S < min_angle) continue; //angle less than 2 mrad
       //if decay products fall into the same strip i.e. dead zones
       //if(abs(x[0][0]-x[1][0]) < same_strip) continue; //300 mum 
       //if(abs(y[0][0]-y[1][0]) < same_strip) continue;

       //minimum distance between tracks, 0 - 29S, 1 - proton should not exceed 180 um
       det01 = (x[0][0]-x[1][0])*(y[0][0]*z[1][0]-y[1][0]*z[0][0]) + (y[0][0]-y[1][0])*(-x[0][0]*z[1][0]+x[1][0]*z[0][0]) + (z[0][0]-z[1][0])*(x[0][0]*y[1][0]-x[1][0]*z[0][0]);
       //det12=(x(1,1)-x(1,2))*(dy(1)*dz(2)-dy(2)*dz(1)) +(y(1,1)-y(1,2))*(-dx(1)*dz(2)+dx(2)*dz(1)) +(z(1,1)-z(1,2))*(dx(1)*dy(2)-dx(2)*dz(1))
       dist01 = abs(det01/(rad[0]*rad[1]*sqrt(1.0-cos_ang*cos_ang)));
       //if(dist01 > max_dist) continue;
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
       //if(s_pz > 35.8755 && s_pz < 35.94) ang_s_p_cut->Fill(ang_p1S);
       ang_s_p->Fill(ang_p1S);
       h_s_pz->Fill(s_pz);
       ang_momz->Fill(ang_p1S,s_pz);
       h_qval2->Fill(q_pcm_p*1000);
       //if(s_pz > 35.8755 && s_pz < 35.94) ang_momz->Fill(s_pz,ang_p1S);
  /*     
       ang_momt->Fill(s_pt,ang_p1S);
       ang_momtp->Fill(p_pt,ang_p1S);
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
*/

       //******************************************************************
       //cout << "Num of prim "<<primary << endl;

       if(traPerEvent>0)         delete[] Tra;
       if(MCtracksPerEvent>0)    delete[] track;

    }//loop for nevents

    TCanvas *c_thes = new TCanvas("thesis", "for thesis",0,900,800,900);
    gStyle->SetOptStat(0);
    c_thes->Divide(2,1);
    c_thes->cd(1);
    h_s_pz->GetXaxis()->CenterTitle();
    h_s_pz->GetXaxis()->SetTitle("p_{||}(HI) (GeV/c)");
    h_s_pz->GetXaxis()->SetLabelSize(0.05);
    h_s_pz->GetXaxis()->SetTitleSize(0.06);
    h_s_pz->GetYaxis()->CenterTitle();
    h_s_pz->GetYaxis()->SetTitle("Intensity (counts)");
    h_s_pz->GetYaxis()->SetLabelSize(0.05);
    h_s_pz->GetYaxis()->SetTitleSize(0.06);
    h_s_pz->Draw();
    c_thes->cd(2);
    ang_s_p->SetFillColor(kBlue-8);
    ang_s_p->GetXaxis()->CenterTitle();
    ang_s_p->GetXaxis()->SetTitle("#theta_{HI-p} (mrad)");
    ang_s_p->GetXaxis()->SetLabelSize(0.05);
    ang_s_p->GetXaxis()->SetTitleSize(0.06);
    ang_s_p->GetYaxis()->CenterTitle();
    ang_s_p->GetYaxis()->SetTitle("Intensity (counts)");
    ang_s_p->GetYaxis()->SetLabelSize(0.05);
    ang_s_p->GetYaxis()->SetTitleSize(0.06);
    ang_s_p->Draw();

    TCanvas *c_thes1 = new TCanvas("thesis1", "for thesis1",0,900,800,900);
    gStyle->SetOptStat(0);
    ang_momz->GetYaxis()->CenterTitle();
    ang_momz->GetYaxis()->SetTitle("p_{||}(HI) (GeV/c)");
    ang_momz->GetYaxis()->SetLabelSize(0.05);
    ang_momz->GetYaxis()->SetTitleSize(0.06);
    ang_momz->GetXaxis()->CenterTitle();
    ang_momz->GetXaxis()->SetTitle("#theta_{HI-p} (mrad)");
    ang_momz->GetXaxis()->SetLabelSize(0.05);
    ang_momz->GetXaxis()->SetTitleSize(0.06);
    ang_momz->GetZaxis()->SetLabelSize(0.05);
    ang_momz->Draw("colz");

    TCanvas *c_thes2 = new TCanvas("thesis2", "for thesis2",0,900,800,900);
    gStyle->SetOptStat(0);
    h_qval2->SetFillColor(kBlue-8);
    h_qval2->GetXaxis()->CenterTitle();
    h_qval2->GetXaxis()->SetTitle("Q (MeV)");
    h_qval2->GetXaxis()->SetLabelSize(0.05);
    h_qval2->GetXaxis()->SetTitleSize(0.06);
    h_qval2->GetYaxis()->CenterTitle();
    h_qval2->GetYaxis()->SetTitle("Intensity (counts)");
    h_qval2->GetYaxis()->SetLabelSize(0.05);
    h_qval2->GetYaxis()->SetTitleSize(0.06);
    h_qval2->Draw();


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

