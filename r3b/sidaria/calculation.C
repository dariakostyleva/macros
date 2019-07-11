#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMath.h"

void calculation(){

	Double_t qvalue, qvalue2;
	Double_t T = 31*0.618; //mother 31Ar
	Double_t mp = 0.938272297;
	Double_t mhi = 29.*0.93149406-0.003156; //30Cl
	Double_t theta = 0.054; //rad

	qvalue = TMath::Sin(theta)*TMath::Sin(theta) * mp*T*(T + 2*(mp + mhi))/(2*mhi*(mp + mhi));
	cout << "qvalue in GeV = "<<qvalue<<endl;

//	Double_t plong = 35.9;
//	qvalue2 = sqrt(plong*plong + mhi*mhi) - mhi - T;
//	cout << "qvalue2 in GeV = "<<qvalue2<<endl;

}