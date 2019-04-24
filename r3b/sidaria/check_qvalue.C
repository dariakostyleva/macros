void check_qvalue()
{

Double_t Pmom = 37.1504;
Double_t Mmom = 27.9509;
Double_t Mp = 0.938272;
Double_t Mhi = 27.0102;
Double_t Pp= 1.31387;
Double_t Phi= 35.8378;



    Double_t qvalue = 0.;



	qvalue = sqrt(Pp*Pp+Mp*Mp) + sqrt(Phi*Phi + Mhi*Mhi) - sqrt(Pmom*Pmom + Mmom*Mmom) - (Mp + Mhi -Mmom);
	cout << "Qvalue = " <<qvalue << endl;


}