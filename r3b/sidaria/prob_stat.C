#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void prob_stat(){
	TFile * f = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_1m.root","READ");
    TH1F * h_s_pz = (TH1F*)f->Get("h_s_pz");

    TFile * f2 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_1k.root","READ");
	TH1F * h_s_pz_2 = (TH1F*)f2->Get("h_s_pz");

	TFile * f3 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_2k.root","READ");
	TH1F * h_s_pz_3 = (TH1F*)f3->Get("h_s_pz");

	TFile * f4 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_700ev.root","READ");
	TH1F * h_s_pz_4 = (TH1F*)f4->Get("h_s_pz");

	TFile * f5 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_100ev.root","READ");
	TH1F * h_s_pz_5 = (TH1F*)f5->Get("h_s_pz");

	TFile * f6 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_300ev.root","READ");
	TH1F * h_s_pz_6 = (TH1F*)f6->Get("h_s_pz");

	TFile * f7 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_200ev.root","READ");
	TH1F * h_s_pz_7 = (TH1F*)f7->Get("h_s_pz");

	TFile * f8 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_150ev.root","READ");
	TH1F * h_s_pz_8 = (TH1F*)f8->Get("h_s_pz");

	TFile * f9 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_400ev.root","READ");
	TH1F * h_s_pz_9 = (TH1F*)f9->Get("h_s_pz");

    TFile * f10 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_250ev.root","READ");
    TH1F * h_s_pz_10 = (TH1F*)f10->Get("h_s_pz");

    TFile * f11 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_180ev.root","READ");
    TH1F * h_s_pz_11 = (TH1F*)f11->Get("h_s_pz");

    TFile * f12 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_220ev.root","READ");
    TH1F * h_s_pz_12 = (TH1F*)f12->Get("h_s_pz");

    TFile * f13 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_210ev.root","READ");
    TH1F * h_s_pz_13 = (TH1F*)f13->Get("h_s_pz");

    TFile * f14 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_10ev.root","READ");
    TH1F * h_s_pz_14 = (TH1F*)f14->Get("h_s_pz");

    TFile * f15 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_80ev.root","READ");
    TH1F * h_s_pz_15 = (TH1F*)f15->Get("h_s_pz");

    TFile * f16 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_20ev.root","READ");
    TH1F * h_s_pz_16 = (TH1F*)f16->Get("h_s_pz");

    TFile * f17 = new TFile("/u/dkostyl/R3BRoot/macros/r3b/sidaria/root_with_hist/out_hist_50ev.root","READ");
    TH1F * h_s_pz_17 = (TH1F*)f17->Get("h_s_pz");

    TCanvas * c1 = new TCanvas();
    c1->Divide(1,8);
    c1->cd(1);
    h_s_pz->Draw();

    c1->cd(2);
    h_s_pz_9->Draw();

    c1->cd(3);
    h_s_pz_6->Draw();

    c1->cd(4);
    h_s_pz_10->Draw();

    c1->cd(5);
    h_s_pz_12->Draw();

    c1->cd(6);
    h_s_pz_7->Draw();

    c1->cd(7);
    h_s_pz_8->Draw();

    c1->cd(8);
    h_s_pz_5->Draw();

    c1->cd(8);
    h_s_pz_15->Draw();

    Double_t ks2 = h_s_pz->KolmogorovTest(h_s_pz_2); 
    Double_t ks3 = h_s_pz->KolmogorovTest(h_s_pz_3); 
    Double_t ks4 = h_s_pz->KolmogorovTest(h_s_pz_4); 
	Double_t ks5 = h_s_pz->KolmogorovTest(h_s_pz_5); 
	Double_t ks6 = h_s_pz->KolmogorovTest(h_s_pz_6); 
	Double_t ks7 = h_s_pz->KolmogorovTest(h_s_pz_7);
	Double_t ks8 = h_s_pz->KolmogorovTest(h_s_pz_8);
	Double_t ks9 = h_s_pz->KolmogorovTest(h_s_pz_9);
    Double_t ks10 = h_s_pz->KolmogorovTest(h_s_pz_10);
    Double_t ks11 = h_s_pz->KolmogorovTest(h_s_pz_11);
    Double_t ks12 = h_s_pz->KolmogorovTest(h_s_pz_12);
    Double_t ks13 = h_s_pz->KolmogorovTest(h_s_pz_13);
    Double_t ks14 = h_s_pz->KolmogorovTest(h_s_pz_14);
    Double_t ks15 = h_s_pz->KolmogorovTest(h_s_pz_15);
    Double_t ks16 = h_s_pz->KolmogorovTest(h_s_pz_16);
    Double_t ks17 = h_s_pz->KolmogorovTest(h_s_pz_17);
    printf("KolmogorovTest = %f , 10e6 vs 2k events\n",ks3);
    printf("KolmogorovTest = %f , 10e6 vs 1k events\n",ks2);
    printf("KolmogorovTest = %f , 10e6 vs 700 events\n",ks4);
    printf("KolmogorovTest = %f , 10e6 vs 400 events\n",ks9);
    printf("KolmogorovTest = %f , 10e6 vs 300 events\n",ks6);
    printf("KolmogorovTest = %f , 10e6 vs 200 events\n",ks7);    
    printf("KolmogorovTest = %f , 10e6 vs 150 events\n",ks8); 
    printf("KolmogorovTest = %f , 10e6 vs 100 events\n",ks5);
    printf("KolmogorovTest = %f , 10e6 vs 250 events\n",ks10);
    printf("KolmogorovTest = %f , 10e6 vs 220 events\n",ks12);
    printf("KolmogorovTest = %f , 10e6 vs 210 events\n",ks13);
    printf("KolmogorovTest = %f , 10e6 vs 80 events\n",ks15);
    printf("KolmogorovTest = %f , 10e6 vs 10 events\n",ks14);
    printf("KolmogorovTest = %f , 10e6 vs 20 events\n",ks16);
    printf("KolmogorovTest = %f , 10e6 vs 50 events\n",ks17);
    printf("events in 1m %d\n",(Int_t)h_s_pz->GetEntries());

    TCanvas *c2 = new TCanvas();
    const int n_points = 10;
    Int_t x_events[n_points]={
                              //(Int_t)h_s_pz_14->GetEntries(),
                              //(Int_t)h_s_pz_16->GetEntries(),
                              //(Int_t)h_s_pz_17->GetEntries(),
                              (Int_t)h_s_pz_15->GetEntries(),
                              (Int_t)h_s_pz_5->GetEntries(), 
                              (Int_t)h_s_pz_8->GetEntries(), 
                              (Int_t)h_s_pz_11->GetEntries(),
                              (Int_t)h_s_pz_7->GetEntries(),
                              (Int_t)h_s_pz_13->GetEntries(),
                              (Int_t)h_s_pz_12->GetEntries(),
                              (Int_t)h_s_pz_10->GetEntries(),
                              (Int_t)h_s_pz_6->GetEntries(), 
                              (Int_t)h_s_pz_9->GetEntries(),
                              //(Int_t)h_s_pz_4->GetEntries()
                            };
    Int_t y_prob[n_points]  ={//(Int_t)(ks14*100.),
                              //(Int_t)(ks16*100.),
                              //(Int_t)(ks17*100.),
                              (Int_t)(ks15*100.),
                              (Int_t)(ks5*100.), 
                              (Int_t)(ks8*100.), 
                              (Int_t)(ks11*100.),
                              (Int_t)(ks7*100.),
                              (Int_t)(ks13*100.),
                              (Int_t)(ks12*100.),
                              (Int_t)(ks10*100.),
                              (Int_t)(ks6*100.), 
                              (Int_t)(ks9*100.),
                              //(Int_t)(ks4*100.),
                              };
    TGraph * graph = new TGraph(n_points,x_events,y_prob);
    graph->SetTitle("Probability dependence on entries number");
    graph->SetMarkerSize(2);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlack);
    graph->SetLineColor(kBlue);
    graph->GetXaxis()->SetTitle("Number of entries");
    graph->GetYaxis()->SetTitle("Probability from KS test (%)");
    graph->Draw("AP");
    TF1 *fit1 = new TF1("fit1","pol1",50,130);
    TF1 *fit2 = new TF1("fit2","pol3",110,230);
    TF1 *fit3 = new TF1("fit3","pol1",220,310);
    
    //graph->Fit("fit1","R");
   // graph->Fit("fit2","R+");
   // graph->Fit("fit3","R+");
    TF1 *total = new TF1("total","fit1+fit2+fit3",50,310);

    Double_t par[8];
    fit1->GetParameters(&par[0]);
    fit2->GetParameters(&par[2]);
    fit3->GetParameters(&par[6]);
    total->SetParameters(par);
    total->SetLineColor(kRed);
  //  graph->Fit("total","R");

    //total = new TF1("total","fit1+fit2",60,300);
    TLine *mid = new TLine(50,50,300,50);
    mid->SetLineWidth(4);
    mid->SetLineColor(kBlue);
    mid->Draw("same");
    
}