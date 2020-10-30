#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void qval_mom_thes(){

  auto c = new TCanvas("mycanvas1","mycanvas1",0,0,1000,900);
  gStyle->SetOptStat(0);
 // c->Divide(1,2);
	const int n_points = 10;
  Double_t err_p[n_points];
  Double_t err_ql[n_points]; //lower error bar
  Double_t err_qh[n_points]; //upper error bar
  for (auto i=0; i<n_points; i++){
    err_p[i] = 0.0;
    err_ql[i] = 0.0;
    err_qh[i] = 0.0;
    if(i==3) {err_p[i]= 0.0;} 
  }

	Double_t q_vals[n_points]={0.001,    0.0015,   0.002,    0.0025,   0.003,    0.005,    0.007,    0.009,    0.011,    0.015 };
	Double_t p_vals[n_points]={0.427395, 0.417477, 0.483253, 0.461929, 0.498699, 0.581570, 0.619301, 0.638355, 0.668360, 0.729160}; // range for 100 events


	//TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,&err_p[0]);
  TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,0);
	graph->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
	//graph->SetTitleSize(1.);
	graph->SetMarkerSize(1.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
  graph->GetXaxis()->SetLimits(0.0,0.017);
  graph->GetYaxis()->SetRangeUser(0.4,1.0);
  graph->GetXaxis()->SetTitle("Q (GeV)");
  graph->GetYaxis()->SetTitle("#Delta p_{||}(HI) (GeV/c)");
  graph->GetXaxis()->SetLabelSize(0.05);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->SetLabelSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->CenterTitle();

  TF1 *xsq = new TF1("xsq","[0]*sqrt(x) + [1]",0.000,0.02);
  graph->Fit("xsq","W","",0.000,0.02);
  graph->Draw("AP");
  	
  /*  
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

  TLine *acc2 = new TLine(0,0.897532,0.017,0.897532);
  acc2->SetLineWidth(4);
  acc2->SetLineColor(kGreen);

  

 // c->cd(1);
  graph->Draw("AP");
//  acc->Draw("same");
  acc2->Draw("same");
  c->Update();

  TLegend *l = new TLegend(0.1,0.7,0.48,0.9);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l->AddEntry(graph,"data points","lp");
  //l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
  l->AddEntry(acc2,"Super-FRS acceptance, #pm 2.5%","l");
  l->AddEntry(xsq,"Sqare root fit","l");
  l->Draw("same");
  

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
   graph_asym->GetXaxis()->SetLimits(0.000,0.004);
   graph_asym->GetYaxis()->SetRangeUser(0.4,0.7);
   graph_asym->Draw("a2");
   graph_asym->Draw("p");
   xsq->Draw("same");

  TLegend *l2 = new TLegend(0.1,0.7,0.48,0.9);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l2->AddEntry(graph_asym,"data points","lp");
  l2->AddEntry(xsq,"square root fit function","lp");
  //l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
  l2->Draw("same");
  */
}







