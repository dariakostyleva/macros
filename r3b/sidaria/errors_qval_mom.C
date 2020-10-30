#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"

void errors_qval_mom(){
    //graph with errors both on delta p and counts for refrence qvalue
    // this is a non-iformative graph, because on x we put count number but error bars are in MeV
    auto * c = new TCanvas();
    const Int_t n_points = 7;
    Double_t counts[n_points] = 	 {2.0,       5.0,      10.0,     15.0,     20.0,     25.0,     30.0 };
    Double_t mean_p_vals[n_points] = {0.0852584, 0.141054, 0.166669, 0.175131, 0.179343, 0.183679, 0.187606};
    Double_t err_p[n_points] = 		 {0.0352077, 0.0343829,0.0196286,0.0159129,0.0162614,0.0101512,0.00818744};
	Double_t err_ql[n_points] =      {0.636122,  0.622609, 0.369617, 0.302543, 0.308891, 0.195863, 0.15876}; //lower error bar
	Double_t err_qh[n_points] =      {0.757511,  0.738376, 0.407346, 0.32734,  0.334787, 0.205954, 0.165325}; //upper error bar

    auto *graph_asym_full = new TGraphAsymmErrors(n_points,counts,mean_p_vals,&err_ql[0],&err_qh[0],&err_p[0],&err_p[0]);
    graph_asym_full->Draw();
}