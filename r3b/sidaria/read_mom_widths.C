#include<iostream>
#include<fstream>

using namespace std;

void read_mom_widths() {

	//********** Reading values from .txt file *******************************
	// Variable declarations
    const double col_numb = 11;
	fstream in;
	vector < vector <double> > values; // 2d array as a vector of vectors
	vector <double> rows(col_numb); // vector to add into 'array' (represents a row)
	//int row_counter = 0; // Row counter

	// Read file
	in.open("mom_widths.txt"); // Open file
	if (in.is_open()) { // If file has correctly opened...
		// Output debug message
		cout << "File correctly opened" << endl;

		// Dynamically store data into array
		while (in.good()) { // ... and while there are no errors,
			for (int i=0; i < col_numb; i++) {
				in >> rows[i]; // fill the row with col elements
				//cout << "rows[" << i << "] = " << rows[i] << endl;
			}
			values.push_back(rows); // add a new row with vector to array
			//row++; // When for is finished goes to a new row
			if (in.eof()) {
        		break;
            }
		}
	}

	for (int i = 0; i < values.size(); i++) {
      for (int j = 0; j < values[i].size(); j++) {
        cout << "values["<<i<<"]["<<j<<"] = " << values[i][j] << endl;
      }
    }
    //All the values are read into 2d array at this point, where first column is number of events
    //the rest corresponds to delta_p_z for particular number of events
    //*************************************************************************

    //************************ Calculating errors and drawing graph**********************
    vector <double> mean_p; // vector to contain mean of delta_p_z for each number of counts
    vector <double> sum;   // vector to contain sum delta_p_z for each number of counts
    vector <double> std_dev_p;
    

    std::fill(mean_p.begin(), mean_p.end(), 0.0);
    std::fill(sum.begin(), sum.end(), 0.0);
    std::fill(std_dev_p.begin(), std_dev_p.end(), 0.0);

    //nested loop over my big array read from file
	for (int i = 0; i < values.size(); i++) {
      Double_t temp1 = 0.0;
      for (int j = 1; j < values[i].size(); j++) {
      	//cout << "values["<<i<<"]["<<j<<"] = " << values[i][j] << endl;
         temp1 = temp1 + values[i][j];
      }
      sum.push_back(temp1);
      mean_p.push_back(temp1/(values[i].size()-1));
      
      cout << "sum["<<i<<"] = "<< sum[i] << endl;
      cout << "mean_p["<<i<<"] = "<< mean_p[i] << endl;
    }



}