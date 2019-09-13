#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TF1.h"

void pz_hi_depend() {

	TCanvas * c = new TCanvas("c","c",0,0,1000,900);
	c->Divide(2,1);
  c->cd(1);
	const int n_points = 4;
	Double_t thick_t[n_points]={2.6, 1.6 , 0.6, 2.600000e-5};
	Double_t p_vals[n_points]={3.588009, 2.306156, 1.111980, 0.202106};
	TGraph * graph = new TGraph(n_points,thick_t,p_vals);
	graph->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of target thickness");
	//graph->SetTitleSize(1.);
	graph->SetMarkerSize(2.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
	graph->GetXaxis()->SetTitle("Thickness of 9Be target [cm]");
  graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI [GeV/c]");
  graph->Draw("AP");

  TF1 *fit = new TF1("fit","[0]*x + [1]",0,2.8);
  graph->Fit("fit","","",0.0,2.8);

  	
  //line to show FRS acceptance
  TLine *acc = new TLine(0,0.358881,2.8,0.358881);
  acc->SetLineWidth(4);
  acc->SetLineColor(kBlue);

  TLine *acc2 = new TLine(0,0.897532,2.8,0.897532);
  acc2->SetLineWidth(4);
  acc2->SetLineColor(kGreen);

  acc->Draw("same");
  acc2->Draw("same");

  TLegend *l = new TLegend(0.1,0.7,0.48,0.9);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l->AddEntry(graph,"data points","lp");
  l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
  l->AddEntry(acc2,"Super-FRS acceptance, mean of long mom +-2.5%","l");
  l->AddEntry(fit,"Linear fit","l");
  l->Draw("same");
    //c->Update();

  c->cd(2);
  const int n_points2 = 3;
  Double_t q_vals[n_points2]={0.001, 0.002, 0.003 };
  Double_t p_vals2[n_points2]={ 3.472836, 3.588009, 3.582561a };
  TGraph * graph2 = new TGraph(n_points2,q_vals,p_vals2);
  graph2->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
  graph2->SetMarkerSize(2.0);
  graph2->SetMarkerStyle(kFullCircle);
  graph2->SetMarkerColor(kBlack);
  graph2->SetLineColor(kBlue);
  graph2->GetXaxis()->SetTitle("Decay energy [GeV]");
  graph2->GetYaxis()->SetTitle("Range of longitudinal momentum of HI [GeV/c]");
  graph2->Draw("AP");

   

}