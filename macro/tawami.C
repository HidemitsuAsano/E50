const double w = 0.7; //kg weight of fiber (0.5 mm fiber * 140 cm * 2400 )
const double wFe = 7.85/1000./1000.0; //kg/mm^3 SS400
const double g = 9.80665;//acceleration of gravity
const double E = 206000;// N/mm^2 Young's modules of SS400
const double L = 2000.0;// mm  
const double a = 700.0;// mm
const double b = 600.0;// mm
const double W = 10;// mm
//const double H = 100;// mm 

Double_t func_w(Double_t *x,Double_t *par);
Double_t func_self(Double_t *x,Double_t *par);
Double_t func_sum(Double_t *x,Double_t *par);

void tawami(){

  TCanvas *c1 = new TCanvas("c1","c1",800,800);

  TF1* f1sum = new TF1("f1sum",func_sum,100,300,0);
  f1sum->GetXaxis()->SetTitle("H [mm]");
  f1sum->GetYaxis()->SetTitle(" [mm]");
  char title[512];
  sprintf(title,"deflexion of SS400 (W = %0.2f mm)",W);
  f1sum->SetTitle(title);
  f1sum->SetLineColor(3);
  f1sum->GetYaxis()->SetRangeUser(0.001,0.1);
  f1sum->Draw();
  TF1* f1w = new TF1("f1w",func_w,100,300,0);
  f1w->SetLineColor(2);
  f1w->Draw("same");
  TF1* f1self = new TF1("f1self",func_self,100,300,0);
  f1self->SetLineColor(4);
  f1self->Draw("same");
  //f2->Draw("");
  TLegend *leg = new TLegend(0.6,0.70,0.90,0.90);
  leg->AddEntry(f1sum,"by total weight");
  leg->AddEntry(f1w,"by fiber's weight");
  leg->AddEntry(f1self,"by SS400's own weight");
  leg->Draw();
  c1->SaveAs("tawami.pdf","PDF");
  c1->SaveAs("tawami.png","PNG");
}


Double_t func_w(Double_t *x,Double_t *par){

  Double_t val = w*g/384.0/E*(8.*TMath::Power(L,3.) - 4.*L*TMath::Power(b,2.) + TMath::Power(b,3.))/ ( W/12.0 *TMath::Power(x[0],3.));

  return val;

}

Double_t func_self(Double_t *x,Double_t *par){

  Double_t val = 5.*wFe*W*x[0]*g*TMath::Power(L,4.)/384.0/E/ ( W/12.0 *TMath::Power(x[0],3.));

  return val;

}

Double_t func_sum(Double_t *x,Double_t *par){

  Double_t valw = w*g/384.0/E*(8.*TMath::Power(L,3.) - 4.*L*TMath::Power(b,2.) + TMath::Power(b,3.))/ ( W/12.0 *TMath::Power(x[0],3.));
  Double_t valself = 5.*wFe*W*x[0]*g*TMath::Power(L,4.)/384.0/E/ ( W/12.0 *TMath::Power(x[0],3.));

  return valw+valself;
}
