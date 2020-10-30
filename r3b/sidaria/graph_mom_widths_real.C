void graph_mom_widths_real(){
	auto c = new TCanvas("mycanvas1","mycanvas1",0,0,1000,900);
	const int n_points = 35;

	Double_t event[n_points]  ={5,        6,        7,        8,         9,        10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,      22,        23,       24,       25,      26,       27,       28,       29,       35,       /*37,       38,       */39,/*       40,  */     45,       49,      54,       60};//, 100000, 500000};
	Double_t mean[n_points]   ={0.243872, 0.27267,  0.298167, 0.302438,  0.317804, 0.317971, 0.336922, 0.349343, 0.345063, 0.362955, 0.357824, 0.372115, 0.373634, 0.380843, 0.378248, 0.382244, 0.395443, 0.40159,  0.391219, 0.395291, 0.395225,0.397192, 0.391055, 0.405497, 0.407241, 0.418494, /*0.418062, 0.418279, */0.426374, /*0.422065,*/ 0.430488, 0.434215,0.432012, 0.440249};//, 0.550201, 0.565655};
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
  //graph->GetXaxis()->SetLimits(0.0,63);
  graph->GetXaxis()->SetLimits(0.0,60);
  graph->GetYaxis()->SetRangeUser(0.0,0.6);


  TLine *lim = new TLine(0,full_p_range,60,full_p_range);
  lim->SetLineWidth(2);
  lim->SetLineColor(kBlue);
  lim->Draw("same");

  TLegend *l = new TLegend(0.4,0.6,0.89,0.89);
  //legend->SetHeader("","C");         // option "C" allows to center the header
  l->AddEntry(graph,"Data points with statistical error bars","lp");
  l->AddEntry(lim,"Maximum possible range of REAL longitudinal momentum","l");
  l->Draw("same");

  //figure to show random error on large scale
  const int n_points_t = 8;
  Double_t event_t[n_points]  =   {5,        10,       35,        60,        100,          600,           1000,        10000};
  Double_t mean_rand[n_points]   ={0.243872, 0.317971, 0.418494,  0.440249,  0.4304408333, 0.50634975,   0.520057,    0.542563};
  Double_t std_dev_rand[n_points]={0.0903748,0.0725565,0.0443915, 0.0352723, 0.0272334944, 0.0166425616, 0.0126615784,0.0108987021};

 
 

  auto c3 = new TCanvas("c3","c3",0,0,1300,1000);
  TGraph * g3 = new TGraph(n_points_t,event_t,std_dev_rand);
  c3->SetLogx();
  g3->SetTitle("Standard deviation vs number of events");
  g3->SetMarkerSize(2);
  g3->SetMarkerStyle(kFullCircle);
  g3->SetMarkerColor(kBlack);
  g3->SetLineWidth(2);
  g3->SetLineColor(kRed);
  g3->GetXaxis()->SetTitle("Number of events");
  g3->GetYaxis()->SetTitle("Standard deviation (GeV/c)");
  g3->GetXaxis()->CenterTitle();
  g3->GetYaxis()->CenterTitle();
  g3->GetXaxis()->SetTitleSize(0.06);
  g3->GetYaxis()->SetTitleSize(0.06);
  g3->GetXaxis()->SetLabelSize(0.05);
  g3->GetYaxis()->SetLabelSize(0.05);
  g3->Draw("apl");

  Double_t mean_sys[n_points]   ={0};
  Double_t std_dev_sys[n_points]={0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017, 0.017};
  //Double_t std_dev_sys[n_points]={0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009};

  TGraph * g4 = new TGraph(n_points_t,event_t,std_dev_sys);
  g4->SetTitle("Standard deviation vs number of events");
  g4->SetMarkerSize(2);
  g4->SetMarkerStyle(kFullCircle);
  g4->SetMarkerColor(kRed);
  g4->Draw("pl");

  TLegend *l2 = new TLegend(0.1,0.7,0.48,0.9);
  l2->AddEntry(g3,"Random error","p");
  //l->AddEntry(acc,"FRS acceptance, mean of long mom +-1%","l");
  l2->AddEntry(g4,"Error reflecting resolution of separator","p");
  l2->Draw("same");



  // making large-scale graph
  //********** Reading values from .txt file *******************************
/*  fstream in; //file we plant to read from
  Int_t col_numb = 2;
  vector < vector <double> > values; // 2d array as a vector of vectors
  vector <double> rows(col_numb); // vector to add into "values" (represents a row with col_number elements)
  vector <double> evt_num;  // number of counts/events
  vector <double> momentum;  // 
  std::fill(evt_num.begin(), evt_num.end(), 0.0);
  std::fill(momentum.begin(), momentum.end(), 0.0);
  // Read file
  in.open("./mom_widths_real_thesis.txt"); // Open file
  if (in.is_open()) { // If file is correctly opened...
    // Output debug message
    cout << "File correctly opened" << endl;

    // Dynamically store data into array
    while (in.good()) { // ... and while there are no errors,
      for (int i=0; i < col_numb; i++) {
        in >> rows[i]; // fill the row with col elements
        //cout << "rows[" << i << "] = " << rows[i] << endl;
      }
      values.push_back(rows); // add a new row with vector to 2d vector
      if (in.eof()) {
            break;
            }
    }
  }
  Int_t dim = 0;        // number of points in graph
  dim = values.size();
  for (int i = 0; i < values.size(); i++) {
      Double_t temp1 = 0.0;
      for (int j = 1; j < values[i].size(); j++) {
         //cout << "values["<<i<<"]["<<j<<"] = " << values[i][j] << endl;
         temp1 = temp1 + values[i][j];
      }
      evt_num.push_back(values[i][0]);
      momentum.push_back(values[i][1]);
      //cout << "evt_num["<<i<<"] = "<< evt_num[i] << endl;
      //cout << "sum["<<i<<"] = "<< sum[i] << endl;
      //cout << "mean_p["<<i<<"] = "<< mean_p[i] << endl;
    }  */

/*
  auto c2 = new TCanvas("c2","c2",0,0,1000,900);
  TGraph * g = new TGraphErrors(dim,&evt_num[0],&momentum[0],0,0);
  g->SetTitle("Range of REAL longitudinal momentum distribution of HI from number of events in simulation");
  g->SetMarkerSize(2);
  g->SetMarkerStyle(kFullCircle);
  g->SetMarkerColor(kBlack);
  g->SetLineWidth(2);
  g->SetLineColor(kRed);
  g->GetXaxis()->SetTitle("Number of events");
  g->GetYaxis()->SetTitle("Range of longitudinal momentum of [GeV/c]");
  g->Draw("ap");

 // g->GetXaxis()->SetLimits(0.0,60);
 // g->GetYaxis()->SetRangeUser(0.0,0.6);
*/


}
