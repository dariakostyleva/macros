#include <iostream>
#include <fstream>
#include <string>
using namespace std;
// stupid script to read file with data directly from simulation BEFORE all the material and tracking
void draw_hist(){
  //string line;
  ifstream myfile ("stupid.txt");
  double a;
  TH1F *h = new TH1F("h","h",100,-0.5,0.5);
  while (myfile>>a)
  {
  	//printf("%f \n", a);
  	h->Fill(a);
  }
  h->Draw();

}