//#include "hist.cc"

void DrawHist(const char *filename="20170615T144153_000015_0.root")
{
  TFile *fin = new TFile(filename,"READ");
  if(!fin) return;
  gStyle->SetOptStat("e");
  //gStyle->SetStatX(0.5);
  //gStyle->SetStatY(0.5);
  const int maxhist = 32;//maxhist
  TH1F *hadcmean_low = new TH1F("hadcmean_low","ADC LOW mean",maxhist,0,maxhist);
  TH1F *hadcmean_high = new TH1F("hadcmean_high","ADC HIGH mean",maxhist,0,maxhist);
  TH1F *hadcstddev_low = new TH1F("hadcstddev_low","ADC LOW std.dev",maxhist,0,maxhist);
  TH1F *hadcstddev_high = new TH1F("hadcstddev_high","ADC HIGH std.dev",maxhist,0,maxhist);
 
  
  /*
  TCanvas *cadclow = new TCanvas("cadclow","ADC low gain");
  TH1I *hadclow;
  cadclow->Divide(8,4);
  for(int ich = 0; ich < maxhist; ich++){
    cadclow->cd(ich+1);
    char histname[256];
    sprintf(histname,"ADC_LOW_%d",ich);
    hadclow = (TH1I*)fin->Get(histname);
    hadclow->GetXaxis()->SetRangeUser(750,850);
    hadclow->Draw();
    double mean_low = hadclow->GetMean();
    hadcmean_low->Fill(ich,mean_low);
    double stddev_low = hadclow->GetRMS();
    hadcstddev_low->Fill(ich,stddev_low);
  }*/

  TCanvas *cadchigh0 = new TCanvas("cadchigh0","ADC high gain0",1000,1000);
  TCanvas *cadchigh1 = new TCanvas("cadchigh1","ADC high gain1",1000,1000);
  TCanvas *cadchigh2 = new TCanvas("cadchigh2","ADC high gain2",1000,1000);
  TCanvas *cadchigh3 = new TCanvas("cadchigh3","ADC high gain3",1000,1000);
  TH1I *hadchigh;
  cadchigh0->Divide(8,4);
  cadchigh1->Divide(8,4);
  cadchigh2->Divide(8,4);
  cadchigh3->Divide(8,4);
  for(int ich = 0; ich < maxhist*4; ich++){
    if(ich< 32) cadchigh0->cd(ich+1);
    else if(32<= ich && ich < 64) cadchigh1->cd(ich+1-32);
    else if(64<= ich && ich < 96) cadchigh2->cd(ich+1-64);
    else if(96<= ich && ich < 128) cadchigh3->cd(ich+1-96);
    
    char histname[256];
    sprintf(histname,"ADC_HIGH_%d",ich);
    hadchigh = (TH1I*)fin->Get(histname);
    hadchigh->GetXaxis()->SetRangeUser(700,1600);
    gPad->SetLogy();
    hadchigh->Draw();
  }

  TCanvas *ctot0 = new TCanvas("ctot0","ctot0",1000,1000);
  TCanvas *ctot1 = new TCanvas("ctot1","ctot1",1000,1000);
  TCanvas *ctot2 = new TCanvas("ctot2","ctot2",1000,1000);
  TCanvas *ctot3 = new TCanvas("ctot3","ctot3",1000,1000);
  TH1I* htot;
  ctot0->Divide(8,4);
  ctot1->Divide(8,4);
  ctot2->Divide(8,4);
  ctot3->Divide(8,4);
  for(int ich=0;ich<maxhist*4;ich++){
    if(ich< 32) ctot0->cd(ich+1);
    else if(32<= ich && ich < 64) ctot1->cd(ich+1-32);
    else if(64<= ich && ich < 96) ctot2->cd(ich+1-64);
    else if(96<= ich && ich < 128) ctot3->cd(ich+1-96);

    char histname[256];
    sprintf(histname,"TDC_TOT_%d",ich);
    htot = (TH1I*)fin->Get(histname);
    htot->GetXaxis()->SetRangeUser(0,200);
    gPad->SetLogy();
    htot->Draw();
  }

  /*
  TCanvas *chitproX = new TCanvas("chitproX","hit profile X");
  TCanvas *chitproU = new TCanvas("chitproU","hit profile U");
  TCanvas *chitproV = new TCanvas("chitproV","hit profile V");
  chitproX->Divide(2,2);
  chitproU->Divide(2,2);
  chitproV->Divide(2,2);
  TH1I *hitpro;
  for(int ilr=0;ilr<4;ilr++){
    char histname[256];
    sprintf(histname,"HitproflayerX%d",ilr);
    hitpro = (TH1I*)fin->Get(histname);
    hitpro->GetXaxis()->SetRangeUser(-6.5,1);
    hitpro->GetXaxis()->SetTitle("mm");
    chitproX->cd(ilr+1);
    gPad->SetLogy(0);
    hitpro->SetMinimum(0);
    hitpro->Draw();
    sprintf(histname,"HitproflayerU%d",ilr);
    hitpro = (TH1I*)fin->Get(histname);
    hitpro->GetXaxis()->SetRangeUser(-6.5,1);
    hitpro->GetXaxis()->SetTitle("mm");
    chitproU->cd(ilr+1);
    gPad->SetLogy(0);
    hitpro->SetMinimum(0);
    hitpro->Draw();
    
    sprintf(histname,"HitproflayerV%d",ilr);
    hitpro = (TH1I*)fin->Get(histname);
    hitpro->GetXaxis()->SetRangeUser(-6.5,1);
    hitpro->GetXaxis()->SetTitle("mm");
    chitproV->cd(ilr+1);
    gPad->SetLogy(0);
    hitpro->SetMinimum(0);
    hitpro->Draw();
  }
  
  
  TCanvas *chitprochX = new TCanvas("chitprochX","hit ch. profile X");
  TCanvas *chitprochU = new TCanvas("chitprochU","hit ch. profile U");
  TCanvas *chitprochV = new TCanvas("chitprochV","hit ch. profile V");
  chitprochX->Divide(2,2);
  chitprochU->Divide(2,2);
  chitprochV->Divide(2,2);
  TH1I *hitproch;
  for(int ilr=0;ilr<4;ilr++){
    char histname[256];
    sprintf(histname,"Hitproflayer_chX%d",ilr);
    hitproch = (TH1I*)fin->Get(histname);
    hitproch->GetXaxis()->SetTitle("fiber");
    chitprochX->cd(ilr+1);
    gPad->SetLogy(0);
    hitproch->SetMinimum(0);
    hitproch->Draw();
    sprintf(histname,"Hitproflayer_chU%d",ilr);
    hitproch = (TH1I*)fin->Get(histname);
    hitproch->GetXaxis()->SetTitle("fiber");
    chitprochU->cd(ilr+1);
    gPad->SetLogy(0);
    hitproch->SetMinimum(0);
    hitproch->Draw();
    sprintf(histname,"Hitproflayer_chV%d",ilr);
    hitproch = (TH1I*)fin->Get(histname);
    hitproch->GetXaxis()->SetTitle("fiber");
    chitprochV->cd(ilr+1);
    gPad->SetLogy(0);
    hitproch->SetMinimum(0);
    hitproch->Draw();
  }


  
  TCanvas *chitcorrX = new TCanvas("cadchighX","X-X correlation",0,0,600,600);
  TCanvas *chitcorrU = new TCanvas("cadchighU","U-U correlation",600,0,600,600);
  TCanvas *chitcorrV = new TCanvas("cadchighV","V-V correlation",1200,0,600,600);
  TH2F *hhitcorr;
  chitcorrX->Divide(2,2);
  chitcorrU->Divide(2,2);
  chitcorrV->Divide(2,2);
  for(int ipat=1; ipat<4 ;ipat++){
    char histname[256];
    sprintf(histname,"HitCorrX0_X%d",ipat);
    hhitcorr = (TH2F*)fin->Get(histname);
    hhitcorr->GetXaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetYaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetXaxis()->SetTitle("mm");
    hhitcorr->GetYaxis()->SetTitle("mm");
    gPad->SetLogy(0);
    chitcorrX->cd(ipat);
    hhitcorr->Draw("");
    
    sprintf(histname,"HitCorrU0_U%d",ipat);
    hhitcorr = (TH2F*)fin->Get(histname);
    hhitcorr->GetXaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetYaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetXaxis()->SetTitle("mm");
    hhitcorr->GetYaxis()->SetTitle("mm");
    gPad->SetLogy(0);
    chitcorrU->cd(ipat);
    hhitcorr->Draw("");
    
    
    sprintf(histname,"HitCorrV0_V%d",ipat);
    hhitcorr = (TH2F*)fin->Get(histname);
    hhitcorr->GetXaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetYaxis()->SetRangeUser(-6.5,1);
    hhitcorr->GetXaxis()->SetTitle("mm");
    hhitcorr->GetYaxis()->SetTitle("mm");
    gPad->SetLogy(0);
    chitcorrV->cd(ipat);
    hhitcorr->Draw("");

  }
  */

  /*
  TCanvas *cadcmean = new TCanvas("cadcmean","ADC mean");
  cadcmean->cd();
  hadcmean_low->SetLineColor(2);
  hadcmean_low->SetTitle("ADC Mean");
  hadcmean_low->SetXTitle("ch#");
  hadcmean_low->SetYTitle("ADC Mean (ch)");
  hadcmean_low->SetMinimum(760);
  hadcmean_low->Draw();
 
  hadcmean_high->SetLineColor(4);
  hadcmean_high->Draw("same");

  TCanvas *cadcstddev = new TCanvas("cadcstddev","ADC std.dev");
  cadcstddev->cd();
  hadcstddev_high->SetLineColor(4);
  hadcstddev_high->SetMinimum(0);
  hadcstddev_high->SetTitle("ADC Std.Dev.");
  hadcstddev_high->SetXTitle("ch#");
  hadcstddev_high->SetYTitle("ADC Std.Dev. (ch)");
  hadcstddev_high->Draw("");
  hadcstddev_low->SetLineColor(2);
  hadcstddev_low->Draw("same");
  */
}





