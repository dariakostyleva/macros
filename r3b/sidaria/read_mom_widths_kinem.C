#include<iostream>
#include<fstream>

using namespace std;

void read_mom_widths_kinem(const Int_t col_numb = 0, TString filename = "") {
	if(col_numb == 0 || filename == "") {cout<<"Warning: 1st argument - column number, 2nd - input file name "<<endl; return;}

	//********** Reading values from .txt file *******************************
	fstream in; //file we plant to read from
	vector < vector <double> > values; // 2d array as a vector of vectors
	vector <double> rows(col_numb); // vector to add into "values" (represents a row with col_number elements)
	// Read file
	in.open(filename); // Open file
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
    vector <double> evt_num;	// number of counts/events
    vector <double> mean_p;		// vector to contain mean of delta_p_z for each number of counts
    vector <double> sum;		// vector to contain sum delta_p_z for each number of counts
    vector <double> std_dev_p;	// standard deviation of delta_p_z for each number of counts
    vector <double> error;		// error 
    Int_t dim = 0;				// number of points in graph
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
    auto c = new TCanvas("c","c",800,10,1200,800);
    auto graph = new TGraphErrors(dim,&evt_num[0],&mean_p[0],0,&error[0]);
    graph->SetTitle("Range of longitudinal momentum distribution of HI from number of events in simulation");
    graph->SetMarkerSize(2);
	  graph->SetMarkerStyle(kFullCircle);
	  graph->SetMarkerColor(kBlack);
	  graph->SetLineWidth(2);
	  graph->SetLineColor(kRed);
	  graph->GetXaxis()->SetTitle("Number of events");
  	graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI with errors [GeV/c]");
  	graph->GetYaxis()->SetRangeUser(0.04,0.21);
  	graph->GetXaxis()->SetRangeUser(0.0,35.5);
  	gStyle->SetEndErrorSize(2);
    graph->Draw("ap");
  	//graph->Fit("pol3","","",0.0,35.5);

  	TLine *lim = new TLine(0,full_p_range,35.5,full_p_range);
  	lim->SetLineWidth(2);
  	lim->SetLineColor(kBlue);
  	lim->Draw("same");

  	TLegend *l = new TLegend(0.4,0.6,0.89,0.89);
    //legend->SetHeader("","C");         // option "C" allows to center the header
    l->AddEntry(graph,"Data points with statistical error bars","lp");
    l->AddEntry(lim,"Maximum possible range of longitudinal momentum","l");
    l->Draw("same");

}