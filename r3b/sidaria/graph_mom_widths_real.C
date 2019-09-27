void graph_mom_widths_real(){
	auto c = new TCanvas("mycanvas1","mycanvas1",0,0,1000,900);
	const int n_points = 35;

	Double_t event[n_points]  ={5,        6,        7,        8,         9,        10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,      22,        23,       24,       25,      26,       27,       28,       29,       35,       /*37,       38,       */39,/*       40,  */     45,       49,      54,       60};
	Double_t mean[n_points]   ={0.243872, 0.27267,  0.298167, 0.302438,  0.317804, 0.317971, 0.336922, 0.349343, 0.345063, 0.362955, 0.357824, 0.372115, 0.373634, 0.380843, 0.378248, 0.382244, 0.395443, 0.40159,  0.391219, 0.395291, 0.395225,0.397192, 0.391055, 0.405497, 0.407241, 0.418494, /*0.418062, 0.418279, */0.426374, /*0.422065,*/ 0.430488, 0.434215,0.432012, 0.440249};
	Double_t std_dev[n_points]={0.0903748,0.0815297,0.0767794,0.0716354, 0.0679143,0.0725565,0.0671749,0.0667718,0.0606996,0.0545669,0.0593452,0.0574959,0.0593588,0.0540521,0.0505336,0.0529056,0.0547932,0.0510677,0.0486471,0.0502443,0.049673,0.0471869,0.0505052,0.0536496,0.0602214,0.0443915,/*0.0358135,0.0416574,*/0.0348059,/*0.0365712,*/0.0336424,0.035436,0.0335517,0.0352723};
	Double_t full_p_range = 0.565838; // maximum possible range of plong momentum for particular conditions

	TGraph * graph = new TGraphErrors(n_points,event,mean,0,std_dev);
	graph->SetTitle("Range of REAL longitudinal momentum distribution of HI from number of events in simulation");
    graph->SetMarkerSize(2);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineWidth(2);
	graph->SetLineColor(kRed);
	graph->GetXaxis()->SetTitle("Number of events");
    graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI with errors [GeV/c]");
    graph->Draw("ap");
    graph->GetXaxis()->SetLimits(0.0,63);
    graph->GetYaxis()->SetRangeUser(0.0,0.6);


    TLine *lim = new TLine(0,full_p_range,63,full_p_range);
  	lim->SetLineWidth(2);
  	lim->SetLineColor(kBlue);
  	lim->Draw("same");

  	TLegend *l = new TLegend(0.4,0.6,0.89,0.89);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l->AddEntry(graph,"Data points with statistical error bars","lp");
    l->AddEntry(lim,"Maximum possible range of REAL longitudinal momentum","l");
    l->Draw("same");


}