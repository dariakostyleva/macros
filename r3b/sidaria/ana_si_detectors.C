
//  ---------------------------------------------------------------------------------
//
//   -- General macro for the analysis of radioactive decays using silicon detectors
//      Author: Jose Luis Rodriguez Sanchez, Daria Kostyleva
//      Last Update: 05/03/19
//      Comments: tracking in-flight decay products (30Cl->29S+p for instanse) using 
//      silicon microstrip AMS detectors
//	         
//	
//  ---------------------------------------------------------------------------------
//
//   Usage: 
//      > root -l ana_si_detectors.C
//	         
//	
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


void ana_si_detectors(){

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
    //TString filename = "sim_out_1m.root";
    //TString filename = "/data.local2/G4_sim_momenta/sim_out_1m.root";
    TString filename = "sim_out.root";
    TFile *file0 = TFile::Open(filename);
    TTree* Tree0 = (TTree*)file0->Get("evt");
    Long64_t nevents = Tree0->GetEntries();
    std::cout<<"Number of entries: "<<nevents<<std::endl;

    TFile * out_hist = new TFile("out_hist.root","RECREATE");

    //HISTOGRAMS DEFINITION
    TH1F* h1_T = new TH1F("h1_T","Primary PDG Code",1000,0,400000);
    TH1F* h2_T = new TH1F("h2_T","Primary Energy (MeV)",2000,0,5000000); //Change this maximum energy
    TH1F* h3_T = new TH1F("h3_T","Primary Theta",200,-1.0,1.0);
    TH1F* h4_T = new TH1F("h4_T","Primary Phi",200,-3.2,3.2);

    //momenta
    TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.76,36.05);
    TH1F* h_p_pz = new TH1F("h_p_pz","Longitudinal momentum of proton (GeV/c)",500,1.1,1.4);
    TH1F* h_s_pt = new TH1F("h_s_pt","Transverse momentum of HI (GeV/c)",500,0.0,0.08);
    TH1F* h_p_pt = new TH1F("h_p_pt","Transverse momentum of proton (GeV/c)",500,0.0,0.08);
    TH1F* h_s_px = new TH1F("h_s_px","Px of HI (GeV/c)",500,-0.1,0.1);
    TH1F* h_p_px = new TH1F("h_p_px","Px of proton (GeV/c)",500,-0.1,0.1);
    TH1F* ang_s_p = new TH1F("ang_s_p","Angle between 29S and proton (mrad)",500,40,55);
    TH1F* ang_s_p_cut = new TH1F("ang_s_p_cut","Angle between 29S and proton (mrad) with selecton on pz s",500,40,55);
    TH2F* ang_momz = new TH2F("ang_momz","Correlation between long mom of HI and angle HI-p",500,35.76,36.05,500,0,55);
    TH2F* ang_momt = new TH2F("ang_momt","Correlation between transv mom of HI and angle HI-p",500,0.0,0.08,500,0,55);
    TH2F* ang_momtp = new TH2F("ang_momtp","Correlation between transv mom of proton and angle HI-p",500,0.0,0.08,500,0,55);
    TH2F* h_pt_pt = new TH2F("h_pt_pt","Correlation between transverse momenta of products",500,0.0,0.08,500,0.0,0.08);
    TH2F* h_pz_pz = new TH2F("h_pz_pz","Correlation between longitudinal momenta of products",500,1.1,1.4,500,35.78,36);
    TH2F* h_pt_pz = new TH2F("h_pt_pz","Correlation between long mom of HI and transv mom of proton",500,0.0,0.08,500,35.78,36);
    TH2F* h_pz_pt = new TH2F("h_pz_pt","Correlation between transv mom of HI and long mom of proton",500,1.1,1.4,500,0.0,0.08);
    TAxis* xaxis = new TAxis();
    TH2F* h_corr = new TH2F("h_corr","Correlation plot reflecting entries",125,35.76,36.05,300,1,3000);
    TH2F* h_pz_pt_s = new TH2F("h_pz_pt_s","pz 29S - x, pt 29S - y",500,35.76,36.05,500,0.0,0.08);
    //energies
    TH1F* h_s_en = new TH1F("h_s_en","Kinetic energy of 29S (GeV)",500,17.5,18.5);
    TH1F* h_p_en = new TH1F("h_p_en","Kinetic energy of proton (GeV)",500,0.0,1.0);
    TH1F* h_s_de = new TH1F("h_s_de","Energy loss of 29S (GeV)",500,0.0,0.001);
    TH1F* h_p_de = new TH1F("h_p_de","Energy loss of proton (GeV)",500,0.0,0.001);

    TH1F* h_qval = new TH1F("h_qval","Q-value (GeV)",500,0.001,0.003);
    TH1F* h_qval1 = new TH1F("h_qval1","Q-value (GeV)",500,0.001,0.003);
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
    Double_t cos_ang, ang_p1S;
    Double_t s_pz, p_pz, p_pt, s_pt, s_px, s_py, p_px, p_py;
    Double_t det01, dist01, de_p, de_s;
    

    //***** Values for cheking qvalue reproducement, taken from asciigenerator *******//
    Double_t qvalue_out = 0.;   // I calculate by the formula
    Double_t qvalue_in = 0.00185; // I set initially
    Double_t Mp = 0.938272297;
    Double_t Mhi = 29.*0.93149406-0.003156;
    Double_t Pmom, Phi, Pp;
    Double_t Mmom = Mhi + Mp + qvalue_in;

    //*************************************
    //for(Int_t i=0;i<100;i++){
    for(Int_t i=0;i<nevents;i++){
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

       //loop in MC mother tracks per event. Now we have 29S and proton here
       for(Int_t h=0;h<MCtracksPerEvent;h++){
         // cout << "MCTrack " <<h<<" in event "<<i<<endl; 
		      if(track[h]->GetMotherId()<0){
			       track[h]->GetMomentum(momentum);	
		         h1_T->Fill(track[h]->GetPdgCode()-1000000000);
             h2_T->Fill(track[h]->GetEnergy()*1000.-track[h]->GetMass()*1000.);
			       primary++;	
             if(track[h]->GetPdgCode()==1000160290) {
              // GetMass here is wrong!
               h_s_en->Fill(track[h]->GetEnergy() - track[h]->GetMass());
              // printf("T kin of Sulphur after detectors?  = %f\n", track[h]->GetEnergy() - 27.010172);
             }
             if(track[h]->GetPdgCode()==2212) {
               h_p_en->Fill(track[h]->GetEnergy() - track[h]->GetMass());
               //cout << "Mass of proton  = " << track[h]->GetMass() << endl;
              // cout << "momentum proton phi = " << momentum.Phi() << endl;
             }            
          }   
        }//for MC

       //loop in tracks per event IN DETECTORS i.e. each mother track is going via this 
       for(Int_t h=0;h<traPerEvent;h++){
           
            //GetDetCopyID() - detector number, 1,2,3 or 4
           //******* 29S *******************
           if(Tra[h]->GetPdi()==1000160290 && Tra[h]->GetDetCopyID()==1){
              Tra[h]->MomentumOut(momentum_s);
              s_pt = sqrt(Tra[h]->GetPxOut()*Tra[h]->GetPxOut() + Tra[h]->GetPyOut()*Tra[h]->GetPyOut());
              s_px = Tra[h]->GetPxOut();
              s_pz = Tra[h]->GetPzOut();
              x[0][0] = Tra[h]->GetXOut();
              y[0][0] = Tra[h]->GetYOut();
              z[0][0] = Tra[h]->GetZOut();
              de_s = Tra[h]->GetEloss();
             // printf("de_s = %f\n",de_s);
            //  cout << "Detector " << Tra[h]->GetDetCopyID() << endl;
           }
           //******* proton ****************
           if(Tra[h]->GetPdi()==2212 && Tra[h]->GetDetCopyID()==1){
              Tra[h]->MomentumOut(momentum_p);
              h3_T->Fill(momentum_p.Theta());
              
              p_pz = Tra[h]->GetPzOut();
              p_px = Tra[h]->GetPxOut();
              p_py = Tra[h]->GetPyOut();
              p_pt = sqrt(Tra[h]->GetPxOut()*Tra[h]->GetPxOut() + Tra[h]->GetPyOut()*Tra[h]->GetPyOut());
              h4_T->Fill(momentum_p.Phi());
              x[1][0] = Tra[h]->GetXOut();
              y[1][0] = Tra[h]->GetYOut();
              z[1][0] = Tra[h]->GetZOut(); 
              de_p = Tra[h]->GetEloss();       
             // printf("de_p = %f\n",de_p);       
            }
              
       }//for tracker

       //*********** Calculation of the angles between tracks *************
       //*********** For only 1st silicon detector now, second index 0 ****

       rad[0] = sqrt(x[0][0]*x[0][0] + y[0][0]*y[0][0] + z[0][0]*z[0][0]);
       rad[1] = sqrt(x[1][0]*x[1][0] + y[1][0]*y[1][0] + z[1][0]*z[1][0]);

       // cos of angles = scalar product/their lengths
       cos_ang = (x[0][0]*x[1][0] + y[0][0]*y[1][0] + z[0][0]*z[1][0])/(rad[0]*rad[1]);
       ang_p1S = acos(cos_ang)*1000.;
       if(ang_p1S<2.0) continue;
       //if decay products fall into the same strip i.e. dead zones
       if(abs(x[0][0]-x[1][0]) < 0.03) continue;
       if(abs(y[0][0]-y[1][0]) < 0.03) continue;

       //minimum distance between tracks, 0 - 29S, 1 - proton should not exceed 180 um
       det01 = (x[0][0]-x[1][0])*(y[0][0]*z[1][0]-y[1][0]*z[0][0]) + (y[0][0]-y[1][0])*(-x[0][0]*z[1][0]+x[1][0]*z[0][0]) + (z[0][0]-z[1][0])*(x[0][0]*y[1][0]-x[1][0]*z[0][0]);
       //det12=(x(1,1)-x(1,2))*(dy(1)*dz(2)-dy(2)*dz(1)) +(y(1,1)-y(1,2))*(-dx(1)*dz(2)+dx(2)*dz(1)) +(z(1,1)-z(1,2))*(dx(1)*dy(2)-dx(2)*dz(1))
       dist01 = abs(det01/(rad[0]*rad[1]*sqrt(1.0-cos_ang*cos_ang)));
       if(dist01 > 0.018) continue;
       //******************************************************************

       //cheking q-value reproducement
       Pmom = sqrt(0.618*30*(0.618*30+ 2*(Mhi + Mp + qvalue_in)));
       Phi = sqrt(s_px*s_px + s_py*s_py + s_pz*s_pz);
       Pp = sqrt(p_px*p_px + p_py*p_py + p_pz*p_pz);

       //qvalue_out = sqrt(Pp*Pp+Mp*Mp) + sqrt(Phi*Phi + Mhi*Mhi) - sqrt(Pmom*Pmom + Mmom*Mmom) - (Mp + Mhi -Mmom);
       qvalue_out = sqrt(Pp*Pp+Mp*Mp) + sqrt(Phi*Phi + Mhi*Mhi) - sqrt(Pmom*Pmom + Mmom*Mmom) - (Mp + Mhi -Mmom);
       //cout << "Qvalue = " <<qvalue_out << endl;


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
       h_qval->Fill(qvalue_out);
       h_qval1->Fill(qvalue_out + de_p + de_s);
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


    //MC TRACK CANVAS
    TCanvas* c1 = new TCanvas("MCTrack","MCTrack",0,0,720,900);
    c1->SetFillColor(0);
    c1->SetFrameFillColor(0);
    c1->cd();
    c1->Divide(2,3);
    c1->cd(1);	h1_T->Draw();
    c1->cd(2);	h2_T->Draw();
    c1->cd(3);	h3_T->Draw();
    TLine* line1 = new TLine(0,0,0,0);
    line1->SetLineStyle(2);
    line1->Draw();
    c1->cd(4);	h4_T->Draw();
    c1->cd(5);
    ang_momtp->Draw();
    c1->cd(6);
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
    c6->Divide(2,2);
    c6->cd(1);
    h_s_en->GetXaxis()->SetTitle("Kin energy of HI (GeV)");
    h_s_en->Draw();
    c6->cd(2);
    h_p_en->GetXaxis()->SetTitle("Kin energy of proton (GeV)");
    h_p_en->Draw();
    c6->cd(3);
    
    c6->cd(4);
    h_p_de->GetXaxis()->SetTitle("Energy loss of proton (GeV)");
    h_p_de->Draw();
    c6->cd(5);
   // h_qval->Draw();
    h_qval1->SetLineColor(2);
  //  h_qval1->Draw("same");


    TCanvas * c7 = new TCanvas("qvla", "Enenrgy loss and qvalue",0,900,800,900);
    c7->Divide(1,2);
    c7->cd(1);
    h_qval->GetXaxis()->SetTitle("Qvalue calculated via kinetic energies after tracking (GeV/c)");
    h_qval->Draw();
    c7->cd(2);
    h_s_de->GetXaxis()->SetTitle("Energy loss of HI (GeV)");
    h_s_de->Draw();
    //c7->cd(3);
    //h_qval1->Draw();



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
    h_p_en->Write();
    h_s_en->Write();
    h_p_de->Write();
    h_s_de->Write();
    h_qval->Write();

   // out_hist->Close();
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

