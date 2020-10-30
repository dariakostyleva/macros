void thesis_numevent_real_kinem(){
	auto thesis = new TCanvas("thesis","thesis",0,0,1000,900);
  thesis->Divide(2,1);

  //************************************************************************************************************************
  //kinematic case
  thesis->cd(1);
    //********** Reading values from .txt file *******************************
  fstream in; //file we plant to read from
  Int_t col_numb = 101;
  vector < vector <double> > values; // 2d array as a vector of vectors
  vector <double> rows(col_numb); // vector to add into "values" (represents a row with col_number elements)
  // Read file
  in.open("mom_widths_2_60_100.txt"); // Open file
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

    //All the values are read into 2d array at this point, where first column is number of events
    //the rest corresponds to delta_p_z for particular number of events
    //*******************************************************************************************

    //************************ Calculating errors and drawing graph******************************
    vector <double> evt_num;  // number of counts/events
    vector <double> mean_p;   // vector to contain mean of delta_p_z for each number of counts
    vector <double> sum;    // vector to contain sum delta_p_z for each number of counts
    vector <double> std_dev_p;  // standard deviation of delta_p_z for each number of counts
    vector <double> error;    // error 
    Int_t dim = 0;        // number of points in graph
    Double_t full_p_range = 0.202251; // maximum possible range of plong momentum for particular conditions
    
    std::fill(evt_num.begin(), evt_num.end(), 0.0);
    std::fill(mean_p.begin(), mean_p.end(), 0.0);
    std::fill(sum.begin(), sum.end(), 0.0);
    std::fill(std_dev_p.begin(), std_dev_p.end(), 0.0);
    std::fill(error.begin(), error.end(), 0.0);

    //nested loop over my big 2d vector created from file
    //calculates mean delta p and number of events and fills them into vectors
  for (int i = 0; i < values.size(); i++) {
      Double_t temp1 = 0.0;
      for (int j = 1; j < values[i].size(); j++) {
         //cout << "values["<<i<<"]["<<j<<"] = " << values[i][j] << endl;
         temp1 = temp1 + values[i][j];
      }
      sum.push_back(temp1);
      mean_p.push_back(temp1/(values[i].size()-1));
      evt_num.push_back(values[i][0]);
      //cout << "evt_num["<<i<<"] = "<< evt_num[i] << endl;
      //cout << "sum["<<i<<"] = "<< sum[i] << endl;
      //cout << "mean_p["<<i<<"] = "<< mean_p[i] << endl;
    }

    //same nested loop to calculate standard deviation and error and fill them into vectors
  for (int i = 0; i < values.size(); i++) {
    Double_t temp2 = 0.0;
      for (int j = 1; j < values[i].size(); j++) {
        temp2 = temp2 + (values[i][j]-mean_p[i])*(values[i][j]-mean_p[i]);
      }
      std_dev_p.push_back(sqrt(temp2/(values[i].size()-2)));
      //error.push_back(std_dev_p[i]/sqrt(values[i].size()-1));   
    error.push_back(std_dev_p[i]); 
    //this if is needed for qval_mom to put error bar on mom there
    if(evt_num[i]==30){
      cout << "std_dev_p["<<i<<"] = "<< std_dev_p[i] << " for count number "<< evt_num[i] <<endl;
      cout << "error["<<i<<"] = "<< error[i] << " for count number "<< evt_num[i] <<endl;
      cout << "mean_p["<<i<<"] = "<< mean_p[i] << " for count number "<< evt_num[i] <<endl;
    }
      //cout << "std_dev_p["<<i<<"] = "<< std_dev_p[i] << endl;
      //cout << "error["<<i<<"] = "<< error[i] << endl;
    }
    dim = values.size();

    //drawing graph!
    auto graph = new TGraphErrors(dim,&evt_num[0],&mean_p[0],0,&error[0]);
    graph->SetTitle("Range of longitudinal momentum distribution of HI from number of events in simulation");
    graph->SetMarkerSize(2.6);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(kBlack);
    graph->SetLineWidth(2);
    graph->SetLineColor(kRed);
    graph->GetXaxis()->SetTitle("Number of events");
    graph->GetYaxis()->SetTitle("#Delta p_{||}(HI) (GeV/c)");
    graph->GetYaxis()->SetRangeUser(0.04,0.21);
    graph->GetXaxis()->SetRangeUser(0.0,60.5);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetXaxis()->SetTitleSize(0.06);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->SetLabelSize(0.05);
    graph->GetYaxis()->SetTitleSize(0.06);
    graph->GetYaxis()->CenterTitle();
    gStyle->SetEndErrorSize(2);
    graph->Draw("ap");
    //graph->Fit("pol3","","",0.0,35.5);

    TLine *lim = new TLine(0,full_p_range,60.5,full_p_range);
    lim->SetLineWidth(5);
    lim->SetLineColor(kBlue);
    lim->Draw("same");

    TLegend *l = new TLegend(0.4,0.6,0.89,0.89);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l->AddEntry(graph,"Data points with stat. error bars","lp");
    l->AddEntry(lim,"Maximum range, #Delta p_{||}(HI) = 0.202 GeV/c","l");
    l->Draw("same");

  //************************************************************************************************************************
  //realistic case
  thesis->cd(2);
	const int n_points = 35;

	Double_t event[n_points]  ={5,        6,        7,        8,         9,        10,       11,       12,       13,       14,       15,       16,       17,       18,       19,       20,       21,      22,        23,       24,       25,      26,       27,       28,       29,       35,       /*37,       38,       */39,/*       40,  */     45,       49,      54,       60};
	Double_t mean[n_points]   ={0.243872, 0.27267,  0.298167, 0.302438,  0.317804, 0.317971, 0.336922, 0.349343, 0.345063, 0.362955, 0.357824, 0.372115, 0.373634, 0.380843, 0.378248, 0.382244, 0.395443, 0.40159,  0.391219, 0.395291, 0.395225,0.397192, 0.391055, 0.405497, 0.407241, 0.418494, /*0.418062, 0.418279, */0.426374, /*0.422065,*/ 0.430488, 0.434215,0.432012, 0.440249};
	Double_t std_dev[n_points]={0.0903748,0.0815297,0.0767794,0.0716354, 0.0679143,0.0725565,0.0671749,0.0667718,0.0606996,0.0545669,0.0593452,0.0574959,0.0593588,0.0540521,0.0505336,0.0529056,0.0547932,0.0510677,0.0486471,0.0502443,0.049673,0.0471869,0.0505052,0.0536496,0.0602214,0.0443915,/*0.0358135,0.0416574,*/0.0348059,/*0.0365712,*/0.0336424,0.035436,0.0335517,0.0352723};
	Double_t full_p_range1 = 0.565838; // maximum possible range of plong momentum for particular conditions

	TGraph * graph1 = new TGraphErrors(n_points,event,mean,0,std_dev);
	graph1->SetTitle("Range of REAL longitudinal momentum distribution of HI from number of events in simulation");
  graph1->SetMarkerSize(2.6);
	graph1->SetMarkerStyle(kFullCircle);
	graph1->SetMarkerColor(kBlack);
	graph1->SetLineWidth(2);
	graph1->SetLineColor(kRed);
	graph1->GetXaxis()->SetTitle("Number of events");
  graph1->GetYaxis()->SetTitle("#Delta p_{||}(HI) (GeV/c)");
  graph1->Draw("ap");
  graph1->GetXaxis()->SetLimits(0.0,63);
  graph1->GetYaxis()->SetRangeUser(0.0,0.6);
  graph1->GetXaxis()->SetLabelSize(0.05);
  graph1->GetXaxis()->SetTitleSize(0.06);
  graph1->GetXaxis()->CenterTitle();
  graph1->GetYaxis()->SetLabelSize(0.05);
  graph1->GetYaxis()->SetTitleSize(0.06);
  graph1->GetYaxis()->CenterTitle();


  TLine *lim1 = new TLine(0,full_p_range1,63,full_p_range1);
  lim1->SetLineWidth(5);
  lim1->SetLineColor(kBlue);
  lim1->Draw("same");

  TLegend *l1 = new TLegend(0.4,0.6,0.89,0.89);
    //legend->SetHeader("","C");         // option "C" allows to center the header
  l1->AddEntry(graph,"Data points with stat. error bars","lp");
  l1->AddEntry(lim,"Maximum range, #Delta p_{||}(HI) = 0.566 GeV/c","l");
  l1->Draw("same");


}