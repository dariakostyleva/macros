//  ---------------------------------------------------------------------------------
//  Macro to obtain the spread of longitudinal momentum of HI from p decay.
//  The structure is the same as in ana_si_detectors.C
//  Author: Daria Kostyleva
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


void ana_width_kinem(Int_t iterator = 0, Int_t countum = 0, Int_t runnumb = 0){

    TString file_in = "sim_out_kinem.root";
    TFile *file0 = TFile::Open(file_in);
    TTree* Tree0 = (TTree*)file0->Get("evt");
    Long64_t nevents = Tree0->GetEntries();
    //std::cout<<"Number of entries: "<<nevents<<std::endl;

    //HISTOGRAMS DEFINITION

    TH1F* h_s_pz = new TH1F("h_s_pz","Longitudinal momentum of HI (GeV/c)",500,35.76,36.05);
  
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

    Int_t traPerEvent=0;
    Int_t MCtracksPerEvent=0;

    Double_t s_pz, delta_s_pz;
    Double_t s_pz_max = 35.9, s_pz_min = 35.9;
    //Double_t widths_arr[runnumb]; // runnumb - this many times this script is run 

    //************* MAIN LOOP OVER EVENTS *************************************
    for(Int_t i=0;i<nevents;i++){
       //if(i%10000 == 0) printf("Event:%i\n",i);

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
	
       // *************** LOOP over tracks per event! ****************************************

       for(Int_t h=0;h<traPerEvent;h++){          
           //******* 29S *******************
           if(Tra[h]->GetPdi()==1000160290 && Tra[h]->GetDetCopyID()==1){
              Tra[h]->MomentumOut(momentum_s);
              s_pz = Tra[h]->GetPzOut();
              if(s_pz > s_pz_max) s_pz_max = s_pz;
              if(s_pz < s_pz_min && s_pz > 0.) s_pz_min = s_pz;
              //printf("pz is %f, pz max is %f, pz min is %f\n",s_pz, s_pz_max,s_pz_min);
           }
       }//for tracker
 
       h_s_pz->Fill(s_pz);
     

       if(traPerEvent>0)         delete[] Tra;

    }//loop for nevents
    //TCanvas *c1 = new TCanvas("long_mom", "Longitudinal momentum of HI",0,0,800,900);
    //h_s_pz->Draw();

    delta_s_pz = s_pz_max-s_pz_min;
    //printf("pz max is %f, pz min is %f\n",s_pz_max,s_pz_min);
    //printf("delta p of HI = %f GeV/c \n",delta_s_pz);


    ofstream myfile;
    myfile.open ("mom_widths.txt",std::ios_base::app);
    if(iterator == 0) myfile << countum << " ";
    myfile << delta_s_pz << " ";
    if(iterator == runnumb-1) myfile << "\n";
    myfile.close();
    //widths_arr[iterator] = delta_s_pz;
    printf("Macro ana_width.C finished succesfully.\n");

  /*  if(iterator==runnumb){
      //printf("We got: \n");
      for (Int_t i=0; i<runnumb+1; i++){
        //printf("Width in %d run = %f\n",i,widths_arr[i]);
      }
    }
  */
}