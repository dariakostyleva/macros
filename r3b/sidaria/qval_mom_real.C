#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void qval_mom_real(){

  auto c = new TCanvas("mycanvas1","mycanvas1",0,0,1000,900);
 // c->Divide(1,2);
	const int n_points = 11;
  Double_t err_p[n_points];
  Double_t err_ql[n_points]; //lower error bar
  Double_t err_qh[n_points]; //upper error bar
  for (auto i=0; i<n_points; i++){
    err_p[i] = 0.0;
    err_ql[i] = 0.0;
    err_qh[i] = 0.0;
    if(i==3) {err_p[i]= 0.0348059;} //0.00818744 GEV/C !!!!!! for 39 counts

  }

	Double_t q_vals[n_points]={0.0005,   0.001,    0.0015,   0.002,    0.0025,   0.003,    0.005,    0.007,    0.009,    0.011,    0.015 };
	Double_t p_vals[n_points]={0.454960, 0.500328, 0.534554, 0.552021, 0.573883, 0.604198, 0.647064, 0.746405, 0.757967, 0.814903, 0.901119};
  // for 0.002 MeV qvalue total range is 0.565838 at ~ 390 000 events

	TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,&err_p[0]);
	graph->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
	//graph->SetTitleSize(1.);
	graph->SetMarkerSize(1.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
	graph->GetXaxis()->SetTitle("Qvalue [GeV]");
  graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI [GeV/c]");
  graph->GetXaxis()->SetLimits(0.0,0.02);
  graph->GetYaxis()->SetRangeUser(0.0,1.0);
  	
  	// sqare root func to fit graph
  TF1 *xsq = new TF1("xsq","[0]*sqrt(x) + [1]",0.000,0.02);
  graph->Fit("xsq","W","",0.000,0.02);
  Double_t p0 = xsq->GetParameter(0); //fit parameter for square function
  Double_t p1 = xsq->GetParameter(1); 
  err_qh[3] = ((p_vals[3] + err_p[3] - p1)/p0)*((p_vals[3] + err_p[3] - p1)/p0) - q_vals[3];
  err_ql[3] = -((p_vals[3] - err_p[3] - p1)/p0)*((p_vals[3] - err_p[3] - p1)/p0) + q_vals[3];
  cout << "err_qh[3] = " << err_qh[3] << " GeV "<< endl;
  cout << "err_ql[3] = " << err_ql[3] << " GeV "<< endl;

  //line to show FRS acceptance
  TLine *acc = new TLine(0,0.358881,6.0,0.358881);
  acc->SetLineWidth(4);
  acc->SetLineColor(kBlue);

  TLine *acc2 = new TLine(0,0.897532,10.0,0.897532);
  acc2->SetLineWidth(4);
  acc2->SetLineColor(kGreen);

  

 // c->cd(1);
  graph->Draw("AP");
  acc->Draw("same");
  acc2->Draw("same");
  c->Update();
/*
    TLegend *l = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l->AddEntry(graph,"data points","lp");
    l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
    l->AddEntry(acc2,"Super-FRS acceptance, mean of long mom +-2.5%","l");
    l->AddEntry(xsq,"Sqare root fit","l");
    l->Draw("same");
  */

   auto c2 = new TCanvas("mycanvas2","mycanvas2",0,0,1000,900);
   auto * graph_asym = new TGraphAsymmErrors(n_points,q_vals,p_vals,&err_ql[0],&err_qh[0],&err_p[0],&err_p[0]);
   graph_asym->SetTitle("Zoomed graph with calculated errors on qvalue");
   graph_asym->GetXaxis()->SetTitle("Qvalue [GeV]");
   graph_asym->GetYaxis()->SetTitle("Range of longitudinal momentum of HI [GeV/c]");
   graph_asym->SetMarkerColor(kBlack);
   graph_asym->SetMarkerSize(2);
   graph_asym->SetFillColor(kGreen);
   graph_asym->SetFillStyle(3001);
   graph_asym->SetMarkerStyle(kFullCircle);
   graph_asym->GetXaxis()->SetLimits(0.001,0.003);
   graph_asym->GetYaxis()->SetRangeUser(0.4,0.7);
   graph_asym->Draw("a2");
   graph_asym->Draw("p");
   xsq->Draw("same");
}







