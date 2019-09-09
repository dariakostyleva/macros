#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void qval_mom(){

  	auto c = new TCanvas("mycanvas1","mycanvas1",0,0,1000,900);
  	c->Divide(1,2);
	const int n_points = 21;
	Double_t err_p[21];
	Double_t err_ql[21]; //lower error bar
	Double_t err_qh[21]; //upper error bar
	for (auto i=0; i<21; i++){
		err_p[i] = 0.0;
		err_ql[i] = 0.0;
		err_qh[i] = 0.0;
		//if(i==3) {err_p[i]= 0.0196286;} //this number is error for 10 counts taken from read_mom_width output
		if(i==3) {err_p[i]= 0.00818744;} //0.00818744 GEV/C !!!!!! for 30 counts

	}

	Double_t q_vals[n_points]={0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.0030, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.012, 0.015, 0.017, 0.019, 0.021, 0.023, 0.025, 0.027};
	Double_t p_vals[n_points]={0.100967, 0.143272, 0.174591, 0.202106, 0.224812, 0.247471, 0.284004, 0.318119, 0.348724, 0.375557, 0.402966, 0.426929, 0.449352, 0.492371, 0.550659, 0.586452, 0.620327, 0.652485,0.685322, 0.712490, 0.740807};
	TGraph * graph = new TGraphErrors(n_points,q_vals,p_vals,0,&err_p[0]);
	graph->SetTitle("Longitudinal momentum range (measured by FRS/SFRS) in dependence of decay energy");
	//graph->SetTitleSize(1.);
	graph->SetMarkerSize(1.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
	graph->GetXaxis()->SetTitle("Qvalue [GeV]");
  	graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI [GeV/c]");
  	
  	// sqare root func to fit graph
  	TF1 *xsq = new TF1("xsq","[0]*sqrt(x)",0,10);
  	graph->Fit("xsq","","",0.0,10.0);
  	Double_t p0 = xsq->GetParameter(0); //fit parameter for square function
  	err_qh[3] = ((p_vals[3] + err_p[3])/p0)*((p_vals[3] + err_p[3])/p0) - q_vals[3];
	err_ql[3] = -((p_vals[3] - err_p[3])/p0)*((p_vals[3] - err_p[3])/p0) + q_vals[3];
	cout << "err_qh[3] = " << err_qh[3] << " GeV "<< endl;
	cout << "err_ql[3] = " << err_ql[3] << " GeV "<< endl;
  	//getting the qvalue error out dependent on momentum systematics error


  	//line to show FRS acceptance
  	TLine *acc = new TLine(0,0.358881,6.0,0.358881);
  	acc->SetLineWidth(4);
  	acc->SetLineColor(kBlue);

  	TLine *acc2 = new TLine(0,0.897532,10.0,0.897532);
  	acc2->SetLineWidth(4);
  	acc2->SetLineColor(kGreen);

    //graph->GetXaxis()->SetLimits(0.0,45.0);
    graph->GetXaxis()->SetLimits(0.0,0.006);
    graph->GetYaxis()->SetRangeUser(0.0,0.40);
    c->cd(1);
    graph->Draw("AP");
    acc->Draw("same");
    acc2->Draw("same");
    c->Update();

    TLegend *l = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l->AddEntry(graph,"data points","lp");
    l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
    l->AddEntry(acc2,"Super-FRS acceptance, mean of long mom +-2.5%","l");
    l->AddEntry(xsq,"Sqare root fit","l");
    l->Draw("same");

    TPaveText *pt = new TPaveText(27.34719,0.4730493,43.86147,0.6732546,"br");
    pt->AddText("Reaction: 30Cl -> 29S + p");
    pt->AddText("Ekin(30Cl) = 18.54 GeV");
    pt->AddText("Reference decay energy = 2 MeV");
    pt->Draw("same");

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
    graph_asym->GetYaxis()->SetRangeUser(0.15,0.25);
    graph_asym->Draw("a2");
    graph_asym->Draw("p");
    xsq->Draw("same");

    TLegend *l2 = new TLegend(0.1,0.7,0.48,0.9);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l2->AddEntry(graph_asym,"data points","lp");
    l2->AddEntry(xsq,"Sqare root fit","l");
    l2->Draw("same");

/*
  	TCanvas* mycanvas = new TCanvas("mycanvas","mycanvas",1300,0, 1000,1000);
	mycanvas->Divide(1,3);
	mycanvas->cd(1);
	graph->Draw("ACP");

  	
  	mycanvas->cd(2);
  	//long mom histo from file and rect fit to it
    TFile * f = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_1m.root","READ");
    TH1F * h_s_pz = (TH1F*)f->Get("h_s_pz");
    TF1 * rect = new TF1("rect", "(x>35.800137 && x<36.000565)*[0] ", 35.6,36.2);
    h_s_pz->Fit("rect");
    Double_t chi2 = rect->GetChisquare();
    cout << "chi2 of rect fit = " << chi2 << endl;
    TH1F * rect_hist = (TH1F*)rect->CreateHistogram();
    cout <<"bin number " << rect_hist->GetNbinsX() << endl;
   // std::cout<<"height = "<<rect->GetParameter(0)<<std::endl;
    Double_t half_height = rect->GetParameter(0)/2;
    printf("FWHM of rect fit = %f\n", half_height);
    TF1 * mid = new TF1("mid","3.871425",30,40); 
	//mid->Draw("same");

    mycanvas->cd(3);
    rect_hist->SetFillColor(kBlue);
    rect_hist->Draw("bar");
*/



/*
 
	TFile * f2 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_1k.root","READ");
	TH1F * h_s_pz_2 = (TH1F*)f2->Get("h_s_pz");

	TFile * f3 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_10k.root","READ");
	TH1F * h_s_pz_3 = (TH1F*)f3->Get("h_s_pz");

	TFile * f4 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_500ev.root","READ");
	TH1F * h_s_pz_4 = (TH1F*)f4->Get("h_s_pz");

	TFile * f5 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_100ev.root","READ");
	TH1F * h_s_pz_5 = (TH1F*)f5->Get("h_s_pz");

	TFile * f6 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_300ev.root","READ");
	TH1F * h_s_pz_6 = (TH1F*)f6->Get("h_s_pz");

	TFile * f7 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_200ev.root","READ");
	TH1F * h_s_pz_7 = (TH1F*)f7->Get("h_s_pz");

	TFile * f8 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_150ev.root","READ");
	TH1F * h_s_pz_8 = (TH1F*)f8->Get("h_s_pz");

	TFile * f9 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/out_hist_400ev.root","READ");
	TH1F * h_s_pz_9 = (TH1F*)f9->Get("h_s_pz");

    TCanvas * c1 = new TCanvas();
    c1->Divide(1,8);
    c1->cd(1);
    h_s_pz->Draw();

    c1->cd(2);
    h_s_pz_2->Draw();

    c1->cd(3);
    h_s_pz_3->Draw();

    c1->cd(4);
    h_s_pz_4->Draw();

    c1->cd(5);
    h_s_pz_5->Draw();

    c1->cd(6);
    h_s_pz_6->Draw();

    c1->cd(7);
    h_s_pz_7->Draw();

    c1->cd(8);
    h_s_pz_8->Draw();

    c1->cd(8);
    h_s_pz_9->Draw();

    Double_t ks2 = h_s_pz->KolmogorovTest(h_s_pz_2); 
    Double_t ks3 = h_s_pz->KolmogorovTest(h_s_pz_3); 
    Double_t ks4 = h_s_pz->KolmogorovTest(h_s_pz_4); 
	Double_t ks5 = h_s_pz->KolmogorovTest(h_s_pz_5); 
	Double_t ks6 = h_s_pz->KolmogorovTest(h_s_pz_6); 
	Double_t ks7 = h_s_pz->KolmogorovTest(h_s_pz_7);
	Double_t ks8 = h_s_pz->KolmogorovTest(h_s_pz_8);
	Double_t ks9 = h_s_pz->KolmogorovTest(h_s_pz_9);

    printf("KolmogorovTest = %f , 10e6 vs 10k events\n",ks3);
    printf("KolmogorovTest = %f , 10e6 vs 1k events\n",ks2);
    printf("KolmogorovTest = %f , 10e6 vs 500 events\n",ks4);
    printf("KolmogorovTest = %f , 10e6 vs 400 events\n",ks9);
    printf("KolmogorovTest = %f , 10e6 vs 300 events\n",ks6);
    printf("KolmogorovTest = %f , 10e6 vs 200 events\n",ks7);    
    printf("KolmogorovTest = %f , 10e6 vs 150 events\n",ks8); 
    printf("KolmogorovTest = %f , 10e6 vs 100 events\n",ks5);

*/

}







