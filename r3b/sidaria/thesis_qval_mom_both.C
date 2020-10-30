#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void thesis_qval_mom_both(){

  //**********************************************************************************************************
  //real
	const int n_points = 11;
  Double_t err_p[n_points];
  Double_t err_ql[n_points]; //lower error bar
  Double_t err_qh[n_points]; //upper error bar
  for (auto i=0; i<n_points; i++){
    err_p[i] = 0.0;
    err_ql[i] = 0.0;
    err_qh[i] = 0.0;
    if(i==3) {err_p[i]= 0.0335517;} // GEV/C !!!!!! for 54 counts

  }

	Double_t q_vals[n_points]={0.0005*1000,   0.001*1000,    0.0015*1000,   0.002*1000,    0.0025*1000,   0.003*1000,    0.005*1000,    0.007*1000,    0.009*1000,    0.011*1000,    0.015*1000 };
	Double_t p_vals[n_points]={0.454960, 0.500328, 0.534554, 0.552021, 0.573883, 0.604198, 0.647064, 0.746405, 0.757967, 0.814903, 0.901119};

  // for 0.002 MeV qvalue total range is 0.565838 at ~ 390 000 events

	//TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,&err_p[0]);
  TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,0);
	graph->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
	//graph->SetTitleSize(1.);
	graph->SetMarkerSize(1.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
  graph->GetXaxis()->SetLimits(0.0,4);
  graph->GetYaxis()->SetRangeUser(0.4,1.0);
  graph->GetXaxis()->SetTitle("Q (GeV)");
  graph->GetYaxis()->SetTitle("#Delta #vec p_{||}(HI) (GeV/c)");
  graph->GetXaxis()->SetLabelSize(0.05);
  graph->GetXaxis()->SetTitleSize(0.07);
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->SetLabelSize(0.05);
  graph->GetYaxis()->SetTitleSize(0.07);
  graph->GetYaxis()->CenterTitle();
 //TLatex text;
 //text.DrawLatexNDC(0.1,0.1, "(a)");
  	// sqare root func to fit graph
  TF1 *xsq = new TF1("xsq","[0]*sqrt(x) + [1]",0.000,4);
  graph->Fit("xsq","W","",0.000,4);
  Double_t p0 = xsq->GetParameter(0); //fit parameter for square function
  Double_t p1 = xsq->GetParameter(1); 
  err_qh[3] = ((p_vals[3] + err_p[3] - p1)/p0)*((p_vals[3] + err_p[3] - p1)/p0) - q_vals[3];
  err_ql[3] = -((p_vals[3] - err_p[3] - p1)/p0)*((p_vals[3] - err_p[3] - p1)/p0) + q_vals[3];
  cout << "err_qh[3] = " << err_qh[3] << " GeV "<< endl;

   auto c2 = new TCanvas("mycanvas2","mycanvas2",0,0,1000,900);
   c2->Divide(2,1);
   c2->cd(2);
   auto * graph_asym = new TGraphAsymmErrors(n_points,q_vals,p_vals,&err_ql[0],&err_qh[0],&err_p[0],&err_p[0]);
   graph_asym->SetTitle("(b) simulation of setup responce");
   graph_asym->GetXaxis()->SetTitle("Q (MeV)");
   graph_asym->GetYaxis()->SetTitle("#Delta #vec p_{||}(HI) (GeV/c)");
   graph_asym->SetMarkerColor(kBlack);
   graph_asym->SetMarkerSize(2.6);
   //graph_asym->SetFillColor(kGreen);
   graph_asym->SetFillStyle(3001);
   graph_asym->SetMarkerStyle(kFullCircle);
   graph_asym->GetXaxis()->SetLimits(0.000,4);
   graph_asym->GetYaxis()->SetRangeUser(0.4,0.7);
   graph_asym->GetXaxis()->SetLabelSize(0.05);
   graph_asym->GetXaxis()->SetTitleSize(0.07);
   graph_asym->GetXaxis()->CenterTitle();
   graph_asym->GetYaxis()->SetLabelSize(0.05);
   graph_asym->GetYaxis()->SetTitleSize(0.07);
   graph_asym->GetYaxis()->CenterTitle();
   graph_asym->Draw("a2");
   graph_asym->Draw("p");
   xsq->Draw("same");

  TLegend *l2 = new TLegend(0.1,0.7,0.48,0.9);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l2->AddEntry(graph_asym,"Data points","lp");
  l2->AddEntry(xsq,"Square root fit","lp");
  //l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
  l2->Draw("same");


  //******************************************************************************************************************************
  //kinematic

  const int n_points1 = 21;
  Double_t err_p1[21];
  Double_t err_ql1[21]; //lower error bar
  Double_t err_qh1[21]; //upper error bar
  for (auto i=0; i<21; i++){
    err_p1[i] = 0.0;
    err_ql1[i] = 0.0;
    err_qh1[i] = 0.0;
    //if(i==3) {err_p[i]= 0.0196286;} //this number is error for 10 counts taken from read_mom_width output
    if(i==3) {err_p1[i]= 0.00818744;} //0.00818744 GEV/C !!!!!! for 30 counts

  }

  Double_t q_vals1[n_points1]={0.0005*1000, 0.001*1000, 0.0015*1000, 0.002*1000, 0.0025*1000, 0.0030*1000, 0.004*1000, 0.005*1000, 0.006*1000, 0.007*1000, 0.008*1000, 0.009*1000, 0.010*1000, 0.012*1000, 0.015*1000, 0.017*1000, 0.019*1000, 0.021*1000, 0.023*1000, 0.025*1000, 0.027*1000};
  Double_t p_vals1[n_points1]={0.100967, 0.143272, 0.174591, 0.202106, 0.224812, 0.247471, 0.284004, 0.318119, 0.348724, 0.375557, 0.402966, 0.426929, 0.449352, 0.492371, 0.550659, 0.586452, 0.620327, 0.652485,0.685322, 0.712490, 0.740807};
  TGraph * graph1 = new TGraphErrors(n_points1,q_vals1,p_vals1,0,0);
  graph1->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
  //graph->SetTitleSize(1.);
  graph1->SetMarkerSize(1.5);
  graph1->SetMarkerStyle(kFullCircle);
  graph1->SetMarkerColor(kBlack);
  graph1->SetLineColor(kBlue);
  graph1->GetXaxis()->SetTitle("Q (GeV)");
  graph1->GetYaxis()->SetTitle("#Delta #vec p_{||}(HI) (GeV/c)");
  graph1->GetXaxis()->SetLabelSize(0.05);
  graph1->GetXaxis()->SetTitleSize(0.07);
  graph1->GetXaxis()->CenterTitle();
  graph1->GetYaxis()->SetLabelSize(0.05);
  graph1->GetYaxis()->SetTitleSize(0.07);
  graph1->GetYaxis()->CenterTitle();
    // sqare root func to fit graph
  TF1 *xsq1 = new TF1("xsq1","[0]*sqrt(x)",0,3);
  graph1->Fit("xsq1","","",0.0,3.0);
  Double_t p01 = xsq1->GetParameter(0); //fit parameter for square function
  err_qh1[3] = ((p_vals1[3] + err_p1[3])/p01)*((p_vals1[3] + err_p1[3])/p01) - q_vals1[3];
  err_ql1[3] = -((p_vals1[3] - err_p1[3])/p01)*((p_vals1[3] - err_p1[3])/p01) + q_vals1[3];
  cout << "err_qh[3] = " << err_qh1[3] << " GeV "<< endl;
  cout << "err_ql[3] = " << err_ql1[3] << " GeV "<< endl;


  c2->cd(1);
  auto * graph_asym1 = new TGraphAsymmErrors(n_points1,q_vals1,p_vals1,&err_ql1[0],&err_qh1[0],&err_p1[0],&err_p1[0]);
  graph_asym1->SetTitle("(a) kinematic calculation");
  graph_asym1->GetXaxis()->SetTitle("Q (MeV)");
  graph_asym1->GetYaxis()->SetTitle("#Delta #vec p_{||}(HI) (GeV/c)");
  graph_asym1->GetXaxis()->SetLabelSize(0.05);
  graph_asym1->GetXaxis()->SetTitleSize(0.07);
  graph_asym1->GetYaxis()->SetLabelSize(0.05);
  graph_asym1->GetYaxis()->SetTitleSize(0.07);
  graph_asym1->GetXaxis()->CenterTitle();
  graph_asym1->GetYaxis()->CenterTitle();
  graph_asym1->SetMarkerColor(kBlack);
  graph_asym1->SetMarkerSize(2.6);
  //graph_asym1->SetFillColor(kGreen);
  graph_asym1->SetFillStyle(3001);
  graph_asym1->SetMarkerStyle(kFullCircle);
  graph_asym1->GetXaxis()->SetLimits(1,2.9);
  graph_asym1->GetYaxis()->SetRangeUser(0.15,0.25);
  graph_asym1->Draw("a2");
  graph_asym1->Draw("p");
  xsq1->Draw("same");

  TLegend *l21 = new TLegend(0.1,0.7,0.48,0.9);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l21->AddEntry(graph_asym1,"Data points","lp");
  l21->AddEntry(xsq1,"Sqare root fit","l");
  l21->Draw("same");
}







