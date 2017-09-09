#include <iostream>
#include <string>
#include <fstream>

void rawhitana(string filename="20170615T144153_000015_0.root")
{
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(0);
  const int maxhist = 32;//maxhist per canvas
  //TH1::SetDefaultSumw2();
  TFile *fin = new TFile(filename.c_str());
  if(!fin->IsOpen()){
    return;
  }
  int ipath = filename.find_last_of("/")+1;
  cout << "path" << ipath << endl;
  int iext  = filename.find_last_of(".");
  string fresult = filename.substr(0,iext);
  //string fresult = filename.substr(ipath,iext-ipath);
  fresult +="_ret.txt";
  ofstream fout;
  fout.open(fresult.c_str());

  TCanvas *cadchigh0 = new TCanvas("cadchigh0","ADC high gain0");
  TCanvas *cadchigh1 = new TCanvas("cadchigh1","ADC high gain1");
  TCanvas *cadchigh2 = new TCanvas("cadchigh2","ADC high gain2");
  TCanvas *cadchigh3 = new TCanvas("cadchigh3","ADC high gain3");
  TH1I *hadchigh=NULL;
  cadchigh0->Divide(8,4);
  cadchigh1->Divide(8,4);
  cadchigh2->Divide(8,4);
  cadchigh3->Divide(8,4);
  TF1 *f1=NULL;// = new TF1("f1","gaus",700,900);
  Double_t x[maxhist*4]={0.0};
  Double_t mean[maxhist*4]={0.0};
  Double_t meanerr[maxhist*4]={0.0};
  Double_t sigma[maxhist*4]={0.0};
  Double_t sigmaerr[maxhist*4]={0.0};
  TGraphErrors *grmean[4]={NULL};
  TGraphErrors *grsigma[4]={NULL};
  for(int i=0;i<4;i++){
    grmean[i] = new TGraphErrors();
    grsigma[i] = new TGraphErrors();
  }

  for(int ich = 0; ich < maxhist*4; ich++){
    if(ich< 32) cadchigh0->cd(ich+1);
    else if(32<= ich && ich < 64) cadchigh1->cd(ich+1-32);
    else if(64<= ich && ich < 96) cadchigh2->cd(ich+1-64);
    else if(96<= ich && ich < 128) cadchigh3->cd(ich+1-96);
    
    char histname[256];
    sprintf(histname,"ADC_HIGH_%d",ich);
    hadchigh = (TH1I*)fin->Get(histname);
    //hadchigh->GetXaxis()->SetRangeUser(750,1200);
    //gPad->SetLogy();
    //hadchigh->Rebin(4);
    hadchigh->GetXaxis()->SetRangeUser(650,950);
    hadchigh->Sumw2();
    hadchigh->Draw("HE");
    gPad->SetLogy();
    
    int maxbin = hadchigh->GetMaximumBin();
    int rms    = hadchigh->GetRMS();
    f1 = new TF1(Form("f%d",ich),"gaus",maxbin-3*rms,maxbin+3*rms);
    f1->SetParameter(1,maxbin);
    f1->SetParameter(2,rms);
    hadchigh->Fit(f1,"0q","",maxbin-2*rms,maxbin+2*rms);

    x[ich] = ich;
    mean[ich] = f1->GetParameter(1);
    meanerr[ich] = f1->GetParError(1);
    sigma[ich] = f1->GetParameter(2);
    sigmaerr[ich] = f1->GetParError(2);
    //fout << ich << " " << mean[ich] << " " << sigma[ich] << endl;
    if(ich< maxhist){
      grmean[0]->SetPoint(ich,x[ich],mean[ich]);
      grmean[0]->SetPointError(ich,0,meanerr[ich]);
      grsigma[0]->SetPoint(ich,x[ich],sigma[ich]);
      grsigma[0]->SetPointError(ich,0,sigmaerr[ich]);
    }
    else if(32<= ich && ich < 64){
      grmean[1]->SetPoint(ich-32,x[ich],mean[ich]);
      grmean[1]->SetPointError(ich-32,0,meanerr[ich]);
      grsigma[1]->SetPoint(ich-32,x[ich],sigma[ich]);
      grsigma[1]->SetPointError(ich-32,0,sigmaerr[ich]);
    }
    else if(64<= ich && ich < 96){
      grmean[2]->SetPoint(ich-64,x[ich],mean[ich]);
      grmean[2]->SetPointError(ich-64,0,meanerr[ich]);
      grsigma[2]->SetPoint(ich-64,x[ich],sigma[ich]);
      grsigma[2]->SetPointError(ich-64,0,sigmaerr[ich]);
    }
    else if(96<= ich && ich < 128){
      grmean[3]->SetPoint(ich-96,x[ich],mean[ich]);
      grmean[3]->SetPointError(ich-96,0,meanerr[ich]);
      grsigma[3]->SetPoint(ich-96,x[ich],sigma[ich]);
      grsigma[3]->SetPointError(ich-96,0,sigmaerr[ich]);
    }
  }
  
  iext  = filename.find_last_of(".");
  string picname = filename.substr(0,iext);
  picname+="_adc0.png";
  cadchigh0->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_adc1.png";
  cadchigh1->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_adc2.png";
  cadchigh2->SaveAs(picname.c_str());

  picname = filename.substr(0,iext);
  picname+="_adc3.png";
  cadchigh3->SaveAs(picname.c_str());

  grmean[0]->SetMarkerColor(2);
  grmean[1]->SetMarkerColor(3);
  grmean[2]->SetMarkerColor(4);
  grmean[3]->SetMarkerColor(5);
  grmean[0]->SetMarkerStyle(21);
  grmean[1]->SetMarkerStyle(21);
  grmean[2]->SetMarkerStyle(21);
  grmean[3]->SetMarkerStyle(21);
  grsigma[0]->SetMarkerColor(2);
  grsigma[1]->SetMarkerColor(3);
  grsigma[2]->SetMarkerColor(4);
  grsigma[3]->SetMarkerColor(5);
  grsigma[0]->SetMarkerStyle(20);
  grsigma[1]->SetMarkerStyle(20);
  grsigma[2]->SetMarkerStyle(20);
  grsigma[3]->SetMarkerStyle(20);

  TCanvas *cfitresult = new TCanvas("cfitresult","fit result");
  cfitresult->Divide(2,1);
  cfitresult->cd(1);
  TMultiGraph *mg1 = new TMultiGraph();
  mg1->SetTitle("pedestal mean;readout ch.; ADC (ch)"); 
  mg1->Add(grmean[0],"p");
  mg1->Add(grmean[1],"p");
  mg1->Add(grmean[2],"p");
  mg1->Add(grmean[3],"p");
  mg1->Draw("a");
  
  cfitresult->cd(2);
  TMultiGraph *mg2 = new TMultiGraph();
  mg2->SetTitle("pedestal sigma;readout ch.; ADC (ch)"); 
  mg2->Add(grsigma[0],"p");
  mg2->Add(grsigma[1],"p");
  mg2->Add(grsigma[2],"p");
  mg2->Add(grsigma[3],"p");
  mg2->Draw("a");
  
  picname = filename.substr(0,iext);
  picname+="_fitret.png";
  cfitresult->SaveAs(picname.c_str());

  
  /*
  //pedestal suppression and determination of hit threshold 
  TCanvas *chit0 = new TCanvas("chit0","ADC high gain0",1500,800);
  TCanvas *chit1 = new TCanvas("chit1","ADC high gain1",1500,800);
  TCanvas *chit2 = new TCanvas("chit2","ADC high gain2",1500,800);
  TCanvas *chit3 = new TCanvas("chit3","ADC high gain3",1500,800);
  //TCanvas *cex = new TCanvas("cex","ADC high gain");
  chit0->Divide(8,4);
  chit1->Divide(8,4);
  chit2->Divide(8,4);
  chit3->Divide(8,4);
  */
  
  TCanvas *chit[8];
  for(int ic=0;ic<8;ic++){
    chit[ic] = new TCanvas(Form("chit%d",ic),"ADC high gain",1500,800);
    chit[ic]->Divide(4,4);
  }

  TGraphErrors *grhitmean[4];
  TGraphErrors *grhitsigma[4];
  for(int i=0;i<4;i++){
    grhitmean[i] = new TGraphErrors();
    grhitsigma[i] = new TGraphErrors();
  }
  
  gStyle->SetOptStat("e");
  hadchigh=NULL;
  TH1I *hhit=NULL;
  Double_t hitmean[maxhist*4];
  Double_t hitmeanerr[maxhist*4];
  Double_t hitsigma[maxhist*4];
  Double_t hitsigmaerr[maxhist*4];
  for(int ich = 0; ich < maxhist*4; ich++){
    int canvasid = (int)ich/16;
    int padid = (int) ich%16;
    chit[canvasid]->cd(1+padid);
    /*
    if(ich< 16) chit->cd(ich+1);
    else if(32<= ich && ich < 64) chit1->cd(ich+1-32);
    else if(64<= ich && ich < 96) chit2->cd(ich+1-64);
    else if(96<= ich && ich < 128) chit3->cd(ich+1-96);
    */

    char histname[256];
    sprintf(histname,"ADC_HIGH_%d",ich);
    hadchigh = (TH1I*)fin->Get(histname);
    gPad->SetLogy();
    hadchigh->Rebin(4);
    hadchigh->GetXaxis()->SetRangeUser(650,3000);
    hadchigh->SetLineColor(1);
    //hadchigh->Draw();
    TH1I* htemp = (TH1I*)hadchigh->Clone();
    htemp->Sumw2();
    htemp->Draw("HE");
    hhit = (TH1I*)hadchigh->Clone();
    hhit->GetXaxis()->SetRangeUser(mean[ich]+10.0*sigma[ich],3000);
    hhit->SetFillColor(7);
    hhit->SetLineColor(7);
    hhit->SetFillStyle(3002);
    hhit->Sumw2();
    hhit->Draw("HEsame");
    x[ich] = ich;
    hitmean[ich] = hhit->GetMean();
    hitmeanerr[ich] = hhit->GetMeanError();
    hitsigma[ich] = hhit->GetRMS();
    hitsigmaerr[ich] = hhit->GetRMSError();
    fout  << ich << " " << mean[ich] << " " 
          << sigma[ich] << " " 
          << hitmean[ich] << " " 
          << hitmeanerr[ich] << " " 
          << hitsigma[ich] << " " 
          << hitsigmaerr[ich] << endl;
    
    if(ich< maxhist){
      grhitmean[0]->SetPoint(ich,x[ich],hitmean[ich]);
      grhitmean[0]->SetPointError(ich,0,hitmeanerr[ich]);
      grhitsigma[0]->SetPoint(ich,x[ich],hitsigma[ich]);
      grhitsigma[0]->SetPointError(ich,0,hitsigmaerr[ich]);
    }
    else if(32<= ich && ich < 64){
      grhitmean[1]->SetPoint(ich-32,x[ich],hitmean[ich]);
      grhitmean[1]->SetPointError(ich-32,0,hitmeanerr[ich]);
      grhitsigma[1]->SetPoint(ich-32,x[ich],hitsigma[ich]);
      grhitsigma[1]->SetPointError(ich-32,0,hitsigmaerr[ich]);
    }
    else if(64<= ich && ich < 96){
      grhitmean[2]->SetPoint(ich-64,x[ich],hitmean[ich]);
      grhitmean[2]->SetPointError(ich-64,0,hitmeanerr[ich]);
      grhitsigma[2]->SetPoint(ich-64,x[ich],hitsigma[ich]);
      grhitsigma[2]->SetPointError(ich-64,0,hitsigmaerr[ich]);
    }
    else if(96<= ich && ich < 128){
      grhitmean[3]->SetPoint(ich-96,x[ich],hitmean[ich]);
      grhitmean[3]->SetPointError(ich-96,0,hitmeanerr[ich]);
      grhitsigma[3]->SetPoint(ich-96,x[ich],hitsigma[ich]);
      grhitsigma[3]->SetPointError(ich-96,0,hitsigmaerr[ich]);
    }
  }
  
  /*
  picname = filename.substr(0,iext);
  picname+="_hit0.png";
  chit0->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_hit1.png";
  chit1->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_hit2.png";
  chit2->SaveAs(picname.c_str());

  picname = filename.substr(0,iext);
  picname+="_hit3.png";
  chit3->SaveAs(picname.c_str());
  */
  for(int ic=0;ic<8;ic++){
    picname = filename.substr(0,iext);
    char suffix[128];
    sprintf(suffix,"_hit%d.png",ic);
    picname += suffix;
    chit[ic]->SaveAs(picname.c_str());
  }
  
  grhitmean[0]->SetMarkerColor(2);
  grhitmean[1]->SetMarkerColor(3);
  grhitmean[2]->SetMarkerColor(4);
  grhitmean[3]->SetMarkerColor(5);
  grhitmean[0]->SetMarkerStyle(21);
  grhitmean[1]->SetMarkerStyle(21);
  grhitmean[2]->SetMarkerStyle(21);
  grhitmean[3]->SetMarkerStyle(21);
  grhitsigma[0]->SetMarkerColor(2);
  grhitsigma[1]->SetMarkerColor(3);
  grhitsigma[2]->SetMarkerColor(4);
  grhitsigma[3]->SetMarkerColor(5);
  grhitsigma[0]->SetMarkerStyle(20);
  grhitsigma[1]->SetMarkerStyle(20);
  grhitsigma[2]->SetMarkerStyle(20);
  grhitsigma[3]->SetMarkerStyle(20);
  
  TCanvas *chits = new TCanvas("chits","hit mean and sigma");
  chits->Divide(2,1);
  chits->cd(1);
  TMultiGraph *mg3 = new TMultiGraph();
  mg3->SetTitle("hit mean;readout ch. ; ADC (ch)"); 
  mg3->Add(grhitmean[0],"p");
  mg3->Add(grhitmean[1],"p");
  mg3->Add(grhitmean[2],"p");
  mg3->Add(grhitmean[3],"p");
  mg3->Draw("a");
  
  chits->cd(2);
  TMultiGraph *mg4 = new TMultiGraph();
  mg4->SetTitle("hit sigma;readout ch. ; ADC (ch)"); 
  mg4->Add(grhitsigma[0],"p");
  mg4->Add(grhitsigma[1],"p");
  mg4->Add(grhitsigma[2],"p");
  mg4->Add(grhitsigma[3],"p");
  mg4->Draw("a");
  picname = filename.substr(0,iext);
  picname+="_pedestalsuppresion.png";
  chits->SaveAs(picname.c_str());
  

  TCanvas *ctot0 = new TCanvas("ctot0","ctot0",1000,1000);
  TCanvas *ctot1 = new TCanvas("ctot1","ctot1",1000,1000);
  TCanvas *ctot2 = new TCanvas("ctot2","ctot2",1000,1000);
  TCanvas *ctot3 = new TCanvas("ctot3","ctot3",1000,1000);
  TH1I* htot=NULL;
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
    //htot->Sumw2();
    htot->Draw("HE");
  }
  
  
  TCanvas *ctotadc0 = new TCanvas("ctotadc0","ctotadc0",1000,800);
  TCanvas *ctotadc1 = new TCanvas("ctotadc1","ctotadc1",1000,800);
  TCanvas *ctotadc2 = new TCanvas("ctotadc2","ctotadc2",1000,800);
  TCanvas *ctotadc3 = new TCanvas("ctotadc3","ctotadc3",1000,800);
  TH2I* htotadc=NULL;
  TH1I* htotadc_padc=NULL;
  TH1I* htotadc_ptot=NULL;
  TH1I* hnhitadc = new TH1I("hnhitadc","nhit adc",128,0,128);
  //hnhitadc->Sumw2(false);
  TH1I* hnhittot = new TH1I("hnhittdc","nhit tdc",128,0,128);
  //hnhittot->Sumw2(false);
  ctotadc0->Divide(8,4);
  ctotadc1->Divide(8,4);
  ctotadc2->Divide(8,4);
  ctotadc3->Divide(8,4);
  for(int ich=0;ich<maxhist*4;ich++){
    if(ich< 32) ctotadc0->cd(ich+1);
    else if(32<= ich && ich < 64) ctotadc1->cd(ich+1-32);
    else if(64<= ich && ich < 96) ctotadc2->cd(ich+1-64);
    else if(96<= ich && ich < 128) ctotadc3->cd(ich+1-96);

    char histname[256];
    sprintf(histname,"AdcToTcorr_%d",ich);
    htotadc = (TH2I*)fin->Get(histname);
    //htotadc->GetXaxis()->SetRangeUser(0,200);
    gPad->SetLogz();
    htotadc->GetXaxis()->SetRangeUser(0,3000);
    htotadc->Draw("colz");
    sprintf(histname,"ADC_HIGH_%d",ich);
    htotadc_padc = (TH1I*) htotadc->ProjectionX();
    float adcthre = mean[ich] + 10.*sigma[ich];
    //cout << adcthre << endl;
    int bin1st = htotadc_padc->FindBin((int)adcthre);
    int binend = htotadc_padc->FindBin(4096);
    //cout << bin1st << "  " << binend << endl;
    int nadchit = htotadc_padc->Integral(bin1st,binend);
    htotadc_ptot = (TH1I*) htotadc->ProjectionY();
    int ntothit = htotadc_ptot->Integral(2,130);
    
    std::cout << "ADC:TOT " << nadchit << "  " << ntothit << std::endl;
    hnhitadc->Fill(ich,nadchit);
    hnhittot->Fill(ich,ntothit);
  }
  picname = filename.substr(0,iext);
  picname+="_adctotcorr0.png";
  ctotadc0->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_adctotcorr1.png";
  ctotadc1->SaveAs(picname.c_str());

  picname = filename.substr(0,iext);
  picname+="_adctotcorr2.png";
  ctotadc2->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_adctotcorr3.png";
  ctotadc3->SaveAs(picname.c_str());


  TCanvas *ccomp0 = new TCanvas("ccomp0","ccomp0",1000,800);
  ccomp0->cd();
  hnhitadc->SetXTitle("ch.#");
  hnhitadc->SetMinimum(0);
  for(int ibin=0;ibin<hnhitadc->GetNbinsX();ibin++){
    int cont = hnhitadc->GetBinContent(ibin+1);
    hnhitadc->SetBinError(ibin+1,sqrt(cont));
    cont = hnhittot->GetBinContent(ibin+1);
    hnhittot->SetBinError(ibin+1,sqrt(cont)); 
  }
  hnhitadc->Draw("H");
  hnhittot->SetLineColor(2);
  hnhittot->Draw("Hsame");
  
  picname = filename.substr(0,iext);
  picname+="_comp.png";
  ccomp0->SaveAs(picname.c_str());
  
  TCanvas *ceff = new TCanvas("ceff","ceff",1000,800);
  ceff->cd();
  TH1D *heff = new TH1D("heff","heff",128,0,128);
  heff->Divide(hnhittot,hnhitadc,1.,1.,"b");
  heff->SetXTitle("ch.#");
  heff->Draw("HE");
  gStyle->SetOptFit(1111);
  heff->Fit("pol0");
  TF1 *fpol = (TF1*)gROOT->GetFunction("pol0");
  Double_t avg = fpol->GetParameter(0);
  ofstream favg;
  favg.open("toteffavg.txt",ios::app);
  cout << avg << endl;
  favg << fresult.c_str() << "  " << avg << endl;
  
  picname = filename.substr(0,iext);
  picname+="_ceff.png";
  ceff->SaveAs(picname.c_str());
  
  
  const char XUVorder[12][2] = {"X","U","V","X","U","V","U","V","X","U","V","X"};
  TCanvas *chitprochX = new TCanvas("chitprochX","hit ch. profile X");
  TCanvas *chitprochU = new TCanvas("chitprochU","hit ch. profile U");
  TCanvas *chitprochV = new TCanvas("chitprochV","hit ch. profile V");
  chitprochX->Divide(2,2);
  chitprochU->Divide(2,2);
  chitprochV->Divide(2,2);
  TH1I *hitproch;
  int cidx =0;
  int cidu =0;
  int cidv =0;
  for(int ilr=0;ilr<12;ilr++){
    char histname[256];
    sprintf(histname,"Hitproflayer_ch%s%d",XUVorder[ilr],ilr);
    hitproch = (TH1I*)fin->Get(histname);
    gPad->SetLogy(0);
    hitproch->SetMinimum(0);
    if(!strcmp(XUVorder[ilr],"X")){
      chitprochX->cd(1+cidx);
      hitproch->Draw();
      cidx++;
    }else if(!strcmp(XUVorder[ilr],"U")){
      chitprochU->cd(1+cidu);
      hitproch->Draw();
      cidu++;
    }else if(!strcmp(XUVorder[ilr],"V")){
      chitprochV->cd(1+cidv);
      hitproch->Draw();
      cidv++;
    }
  }
  
  picname = filename.substr(0,iext);
  picname+="_chipprofchX.png";
  chitprochX->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_chipprofchU.png";
  chitprochU->SaveAs(picname.c_str());
  
  picname = filename.substr(0,iext);
  picname+="_chipprofchV.png";
  chitprochV->SaveAs(picname.c_str());
  
  TCanvas *crawhitcorr[3];
  TH2I* hrawhitcorr=NULL;
  for(int ifirst=0;ifirst<3;ifirst++){
    crawhitcorr[ifirst] = new TCanvas(Form("crawhitcorr%d",ifirst),Form("crawhitcorr%d",ifirst),1000,800);
    crawhitcorr[ifirst]->Divide(4,3);
    for(int ilr=0;ilr<12;ilr++){
      char histname[256];
      sprintf(histname,"rawhitcorr%d_%d",ilr,ifirst);
      hrawhitcorr = (TH2I*)fin->Get(histname);
      crawhitcorr[ifirst]->cd(ilr+1);
      hrawhitcorr->Draw("colz");
    }
  }

  picname = filename.substr(0,iext);
  picname +="_rawhitcorrX.png";
  crawhitcorr[0]->SaveAs(picname.c_str());

  picname = filename.substr(0,iext);
  picname +="_rawhitcorrU.png";
  crawhitcorr[1]->SaveAs(picname.c_str());

  picname = filename.substr(0,iext);
  picname +="_rawhitcorrV.png";
  crawhitcorr[2]->SaveAs(picname.c_str());
}



