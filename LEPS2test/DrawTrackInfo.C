#include <iostream>
#include <string>
#include <fstream>


void DrawTrackInfo(string filename="19.root"){

  TFile *fin = new TFile(filename.c_str(),"READ");
  if(!fin) return;
  gStyle->SetOptStat("emr");
  gStyle->SetStatY(0.9);
  
  const char XUVorder[12][2] = {"X","U","V","X","U","V","U","V","X","U","V","X"};
  TCanvas *cresi = new TCanvas("cresi","cresi",1400,800);
  cresi->Divide(6,2);
  TH1F *hresi[12]={NULL};
  for(int ilr=0;ilr<12;ilr++){
    cresi->cd(ilr+1);
    char hname[256];
    sprintf(hname,"hresi%s%d",XUVorder[ilr],ilr);
    hresi[ilr] = dynamic_cast<TH1F*>(fin->Get(hname));
    hresi[ilr]->Draw();
  }
  int iext  = filename.find_last_of(".");
  string picname = filename.substr(0,iext);
  picname +="_resi.png";
  cresi->SaveAs(picname.c_str());

  TCanvas *cchi2 = new TCanvas("cchi2","cchi2",600,600);
  cchi2->cd();
  TH1F *hchi2=NULL;
  hchi2 = (TH1F*)fin->Get("hchi2");
  hchi2->Draw();
  
  iext  = filename.find_last_of(".");
  picname = filename.substr(0,iext);
  picname +="_chi2.png";
  cchi2->SaveAs(picname.c_str());

  TCanvas *cxy0 = new TCanvas("cxy0","cxy0",800,800);
  cxy0->Divide(2,2);
  TH2F *hxy0 = NULL;
  hxy0 = (TH2F*)fin->Get("hxy0");
  cxy0->cd(1);
  hxy0->Draw("colz");
  cxy0->cd(2);
  TH1F* hx0 = (TH1F*) hxy0->ProjectionX();
  hx0->SetTitle("X0 pos.");
  hx0->Draw();
  cxy0->cd(3);
  TH1F* hy0 = (TH1F*) hxy0->ProjectionY();
  hy0->SetTitle("Y0 pos.");
  hy0->Draw();

  picname = filename.substr(0,iext);
  picname +="_xy0.png";
  cxy0->SaveAs(picname.c_str());

  TCanvas *cnhit = new TCanvas("cnhit","cnhit",800,800);
  cnhit->cd();
  TH1I *hnassociate = (TH1I*)fin->Get("hnassociate");
  hnassociate->Draw();

  picname = filename.substr(0,iext);
  picname +="_nssociate.png";
  cnhit->SaveAs(picname.c_str());

  

}
