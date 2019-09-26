#include<iostream>
#include<fstream>

using namespace std;

void read_mom_widths_real(TString filename = "", Int_t event = 0) {
	if(filename == "") {cout<<"Warning: 1st argument - input file name, 2nd - event number"<<endl; return;}

	//********** Reading values from .txt file *******************************
	fstream in; //file we plant to read from
	vector < vector <double> > values; // 2d array as a vector of vectors
	vector <double> rows(2); // vector to add into "values" (represents a row with 2 elements)
	// Read file
	in.open(filename); // Open file
	if (in.is_open()) { // If file is correctly opened...
		// Output debug message
		cout << "File correctly opened for " <<event<<" events" << endl;

		// Dynamically store data into array
		while (in.good()) { // ... and while there are no errors,
			for (int i=0; i < 2; i++) {
				in >> rows[i]; // fill the row with col elements
				//cout << "rows[" << i << "] = " << rows[i] << endl;
			}
			values.push_back(rows); // add a new row with vector to 2d vector
			if (in.eof()) {
        		break;
            }
		}
	}
	
  //All the values are read into 2d array at this point, where first column is number of events,
  //second column - range of pz
  //*******************************************************************************************

  //************************ Calculating mean and error******************************
    
  
  vector <double> event_array;
  std::fill(event_array.begin(), event_array.end(), 0.0);
  int freq = 0; 
  double sum_p = 0.;
  double mean_p = 0.;
  double std_dev_p = 0.;
  double temp = 0.;

  // here I fill the 2d array for particular event number
	for (int i = 0; i < values.size(); i++) {
    for (int j = 1; j < values[i].size(); j++) {
	    if(values[i][0] == event) {
        event_array.push_back(values[i][j]);
        freq++;
      }
    }
  }

  //if this event is present more than 100 times, I calculate std_dev for it
  if(freq > 100) {
    for(int i = 0; i< 100; i++){
      sum_p = sum_p + event_array[i];
    }
    mean_p = sum_p/100.;
    cout << "mean_p = " << mean_p << endl;
    for(int i = 0; i< 100; i++){
      //cout << event_array[i] << endl;
      temp = temp + (event_array[i] - mean_p)*(event_array[i] - mean_p);
    }
    std_dev_p = sqrt(temp/(100-1));
    cout<< "std_dev_p = " << std_dev_p << endl;
  }
  else cout << event << " events are measured less than 100 times!" << endl;




  //cout << "values.size() " << values.size() << endl;
  //cout << "freq of 20 " << freq << endl;
   /*   sum.push_back(temp1);
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
    graph->SetTitle("Range of longitudinal momentum distribution of HI from number of events in REALISTIC simulation.");
    graph->SetMarkerSize(2);
	  graph->SetMarkerStyle(kFullCircle);
	  graph->SetMarkerColor(kBlack);
	  graph->SetLineWidth(2);
	  graph->SetLineColor(kRed);
	  graph->GetXaxis()->SetTitle("Number of events");
  	graph->GetYaxis()->SetTitle("Range of longitudinal momentum of HI with errors [GeV/c]");
  	//graph->GetYaxis()->SetRangeUser(0.04,0.21);
  	//graph->GetXaxis()->SetRangeUser(0.0,35.5);
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
*/
}
