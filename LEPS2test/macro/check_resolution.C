void check_resolution(const char* filename)
{
  gStyle->SetOptStat("emr");
  TFile::Open(filename);
  TCanvas *cdx0 = new TCanvas("cdx0","cdx0");
  tree->Draw("x0diff>>hdx0(100,-1,1)");
  hdx0->SetXTitle("reconstructed x0 - generated x0 (mm)");
  hdx0->Fit("gaus","","",-0.5,0.5);
  

  TCanvas *cdy0 = new TCanvas("cdy0","cdy0");
  tree->Draw("y0diff>>hdy0(100,-1,1)");
  hdy0->SetXTitle("reconstructed y0 - generated y0 (mm)");
  hdy0->Fit("gaus","","",-0.5,0.5);

  TCanvas *craw = new TCanvas("craw","craw");
  craw->Divide(2,1);
  craw->cd(1);
  tree->Draw("x0>>hx0(100,-6,6)");
  tree->Draw("priposx>>hprix(100,-6,6)");
  hprix->SetLineColor(2);
  hx0->SetXTitle("x position at (z=0) (mm)"); 
  hx0->Draw();
  hprix->Draw("same");
  
  craw->cd(2);
  tree->Draw("y0>>hy0(100,-7,7)");
  tree->Draw("priposy>>hpriy(100,-7,7)");
  hpriy->SetLineColor(2);
  hy0->SetXTitle("y position at (z=0) (mm)"); 
  hy0->Draw();
  hpriy->Draw("same");

  TCanvas *cres = new TCanvas("cres","cres");
  cres->cd();
  tree->Draw("res");

  //TCanvas *c_angleres = new TCanvas("c_angleres","c_angleres")
  //
  

  
}
