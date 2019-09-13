void energy_vs_thickness(){
	TCanvas * c = new TCanvas("c","c");
	c->Divide(2,1);
	const int n_points = 2;
	Double_t vz[n_points]={-1.3, 1.3};
	Double_t en[n_points]={0.618*30, 0.519*30};
	TGraph * graph = new TGraph(n_points,vz,en);
    graph->SetMarkerSize(2.0);
	graph->SetMarkerStyle(kFullCircle);
	graph->SetMarkerColor(kBlack);
	graph->SetLineColor(kBlue);
	c->cd(1);
    graph->Draw("AP");
    TF1 *fit = new TF1("fit","[1]*x + [0]",-1.3,1.3);
    graph->Fit("fit","","",-1.4,1.4);

    const int n = 1000;
	Double_t z[n], y[n];

	for(int i=0; i<n; i++){
		z[i] = 2.6*(gRandom->Uniform(0,1) - 0.5); 
		y[i] = -1.14231*z[i] + 17.055; 
		cout << "z "<< i << " " << z[i] <<" , y = "<< y[i]<< endl;
	}
	c->cd(2);
	TGraph * graph2 = new TGraph(n,z,y);
	graph2->SetMarkerSize(1.0);
	graph2->SetMarkerStyle(kFullCircle);
	graph2->SetMarkerColor(kBlue);
//	graph2->SetLineColor(kBlue);
	graph2->Draw("AP");

}