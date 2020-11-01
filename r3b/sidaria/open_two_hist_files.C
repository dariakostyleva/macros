void open_two_hist_files(){

    //illustration of the influence of momentum spread of 2% (3% in the spread of kinetic energy of mother nucleus)

	TString file_k = "../root_with_hist_kinem/out_hist_100k_noloss.root";
    TFile *filek = TFile::Open(file_k);
    TH1F* hist_ang_k = (TH1F*)filek->Get("ang_s_p");
    TH1F* hist_mom_k = (TH1F*)filek->Get("h_s_pz");


    TString file_ks = "../root_with_hist_kinem/out_hist_100k_noloss_spread.root";
    TFile *fileks = TFile::Open(file_ks);
    TH1F* hist_ang_ks = (TH1F*)fileks->Get("ang_s_p");
    TH1F* hist_mom_ks = (TH1F*)fileks->Get("h_s_pz");

	TString file_r = "../root_with_hist_real/out_hist_100k_nospread.root";
    TFile *filer = TFile::Open(file_r);
    TH1F* hist_ang_r = (TH1F*)filer->Get("ang_s_p");
    TH1F* hist_mom_r = (TH1F*)filer->Get("h_s_pz");


    TString file_rs = "../root_with_hist_real/out_hist_100k_spread.root";
    TFile *filers = TFile::Open(file_rs);
    TH1F* hist_ang_rs = (TH1F*)filers->Get("ang_s_p");
    TH1F* hist_mom_rs = (TH1F*)filers->Get("h_s_pz");


	TCanvas *c = new TCanvas();
	c->Divide(2,2);
    c->cd(1);
    hist_ang_k->GetXaxis()->SetTitle("Relative angle (mrad)");
    hist_ang_k->Draw();

    c->cd(2);
    hist_ang_ks->GetXaxis()->SetTitle("Relative angle (mrad)");
    hist_ang_ks->Draw();

	c->cd(3);
    hist_ang_r->GetXaxis()->SetTitle("Relative angle (mrad)");
    hist_ang_r->Draw();

    c->cd(4);
    //hist_ang_rs->SetLineColor(kRed);
    hist_ang_rs->GetXaxis()->SetTitle("Relative angle (mrad)");
    hist_ang_rs->Draw();
    //hist_ang_rs->SetLineColor(kRed);

    TLegend *l = new TLegend(0.1,0.7,0.48,0.9);
    l->AddEntry("hist_ang_r","No mom spread included","l");
    l->AddEntry("hist_ang_rs","With mom spread included","l");
    //l->Draw("same");


    TCanvas *c2 = new TCanvas();
	c2->Divide(2,2);
	c2->cd(1);
	hist_mom_k->GetXaxis()->SetTitle("Long mom of HI, kinem and no mom spread (GeV/c)");
	hist_mom_k->Draw();

	c2->cd(2);
	hist_mom_ks->GetXaxis()->SetTitle("Long mom of HI, kimen and mom spread (GeV/c)");
	hist_mom_ks->Draw();

    c2->cd(3);
    hist_mom_r->GetXaxis()->SetTitle("Long mom of HI, no mom spread (GeV/c)");
    hist_mom_r->Draw();

    c2->cd(4);
    hist_mom_rs->GetXaxis()->SetTitle("Long mom of HI, with mom spread (GeV/c)");
    hist_mom_rs->Draw();

    TCanvas *c3 = new TCanvas();
    c3->Divide(1,2);
    c3->cd(1);
    hist_ang_k->SetLineColor(kRed);
    hist_ang_k->Draw();
    hist_ang_ks->SetLineColor(kBlue);
    hist_ang_ks->Draw("same");

    c3->cd(2);
    hist_ang_r->Draw();
    hist_ang_rs->SetLineColor(kRed);
    hist_ang_rs->Draw("same");

    TCanvas *c4 = new TCanvas();

    hist_mom_r->SetLineColor(kRed);
    hist_mom_r->Draw();
    hist_mom_rs->SetLineColor(kBlue);
    hist_mom_rs->Draw("same");
    

}