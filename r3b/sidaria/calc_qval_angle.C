#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

void calc_qval_angle(){

	Double_t qvalue;
	Double_t T = 30*0.618; //mother 30Cl
	Double_t mp = 0.938272297;
	Double_t mhi = 29.*0.93149406-0.003156; //30Cl
	Double_t theta = 0.0; //rad

	TFile *f = TFile::Open("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_500k_180719.root");
	TH1F* hist = (TH1F*)f->Get("ang_s_p");
	auto * c = new TCanvas();
	hist->Draw();
	hist->GetNbinsX();
	cout << hist->GetNbinsX() << endl;
	Double_t bmax=0.0; 
	Int_t imax = 0.; 
	for (auto i=1; i<hist->GetNbinsX(); i++){
		if(hist->GetBinContent(i) > bmax ) {
			bmax = hist->GetBinContent(i);
			imax = i;
		}
		cout <<i<<" "<< hist->GetBinContent(i) << endl;
	}
	cout << "bmax = " << bmax << " for imax = "<<imax<<endl;
	theta = (0.004*imax+47.0)/1000; //recalculate bin number into radian
	cout << "theta = " << theta << endl;


	qvalue = TMath::Sin(theta)*TMath::Sin(theta) * mp*T*(T + 2*(mp + mhi))/(2*mhi*(mp + mhi));
	cout << "qvalue in GeV = "<<qvalue<<endl;

}