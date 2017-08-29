//QA_gensft
//quick check for BSFT in GEANT4 info. and generate pdf files
#include <iostream>
#include <math.h>
#include <map>
#include <cstring>

#include <TRint.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TSystem.h>

void QA_gensft(const char* filename)
{
  gStyle->SetOptStat("emr");
  TFile *tf = new TFile(filename, "READ");
  TTree *tree = (TTree*)(tf->Get("tree"));
  
  //number of hits/event
  int  bsftnhits;
  //hit layer /hit segment (=fiber/MPPC)
  std::vector<int> *bsftlayer=0;
  std::vector<int> *bsftseg=0;
  //PID on SFT, track No. on SFT
  std::vector<int> *bsftpid=0;
  std::vector<int> *bsftid=0;
  //E deposit on SFT and time ?
  std::vector<double> *bsftedep=0;
  std::vector<double> *bsfttime=0;
  //absolute value of Momentum  
  std::vector<double> *bsftp=0;
  //theta and phi of momentum vector of beam
  std::vector<double> *bstdtheta=0;
  std::vector<double> *bsftphi=0;
  //local x and y position of sft 
  //local x : phi direction in G4Tubs
  //local y : radial direction in G4Tubs
  std::vector<double> *bsftposlx=0;
  std::vector<double> *bsftposly=0;
  //global x, y and z position of sft
  std::vector<double> *bsftposgx=0;
  std::vector<double> *bsftposgy=0;
  std::vector<double> *bsftposgz=0;
  //pathlength , beta of particle
  //std::vector<double> *bsftpath, *bsftbeta;
  //mass of particle
  //std::vector<double> *bsftmass;
  
 
  tree->SetBranchAddress("sftnhits",&bsftnhits);
  tree->SetBranchAddress("sftlayer",&bsftlayer);
  tree->SetBranchAddress("sftseg",&bsftseg);
  tree->SetBranchAddress("sftp",&bsftp);
  tree->SetBranchAddress("sftposlx",&bsftposlx);
  tree->SetBranchAddress("sftposly",&bsftposly);
  tree->SetBranchAddress("sftposgx",&bsftposgx);
  tree->SetBranchAddress("sftposgy",&bsftposgy);
  tree->SetBranchAddress("sftposgz",&bsftposgz);
  
  TH1I *hsftnhits = new TH1I("hsftnhits","sft nhits",30,0,30);
  TH1I *hsftlayer = new TH1I("hsftlayer","sft hit layer",12,0,12);
  TH1I *hsftseg   = new TH1I("hsftseg","sft hit ch",200,0,200);
  TH1D *hsftp     = new TH1D("hsftp","beam momentum at SFT" ,200,19.5,20.5);
  TH1D *hsftposlx = new TH1D("hsftposlx" ,"local position x ",200,-0.6,0.6);
  TH1D *hsftposly = new TH1D("hsftposly" ,"local position y ",100,-0.1,0.5);   
  TH1D *hsftposgx = new TH1D("hsftposgx", "global position x",600,-60,60); 
  TH1D *hsftposgy = new TH1D("hsftposgy", "global position y",600,-60,60); 
  TH1D *hsftposgz = new TH1D("hsftposgz", "global position z",400,0,400);
  

  Long64_t nentries = tree->GetEntries();
  Long64_t nbytes = 0;
  std::cout << "# of entries " << nentries << std::endl;
  for(Long64_t i=0;i < nentries;i++){
    nbytes += tree->GetEntry(i);
    hsftnhits->Fill(bsftnhits);
    if(bsftnhits != (int)bsftlayer->size()){
      std::cout << "something is wrong" << std::endl;
    }
    if(bsftnhits != (int)bsftseg->size()){
      std::cout << "something is wrong" << std::endl;
    }
    if(bsftnhits != (int)bsftp->size()){
      std::cout << "something is wrong" << std::endl;
    }


    for(int ihit=0; ihit< bsftnhits; ++ihit){
      int sfthitlayer = bsftlayer->at(ihit);
      hsftlayer->Fill(sfthitlayer);
      int sfthitsegment = bsftseg->at(ihit);
      hsftseg->Fill(sfthitsegment);
      double p = bsftp->at(ihit);
      hsftp->Fill(p);
      double poslx = bsftposlx->at(ihit);
      hsftposlx->Fill(poslx);
      double posly = bsftposly->at(ihit);
      hsftposly->Fill(posly);
      double posgx = bsftposgx->at(ihit);
      hsftposgx->Fill(posgx);
      double posgy = bsftposgy->at(ihit);
      hsftposgy->Fill(posgy);
      double posgz = bsftposgz->at(ihit);
      hsftposgz->Fill(posgz);
    }
     
    // hsftlayer->Fill(&bsftlayer); 
  }

  TCanvas *csftnhits = new TCanvas("csftnhits","sft nhits",0,0,400,400);
  csftnhits->cd();
  hsftnhits->SetXTitle("# of nhits on sft");
  hsftnhits->Draw();

  TCanvas *csftlayer = new TCanvas("csftlayer","sft layer",400,0,400,400);
  csftlayer->cd();
  hsftlayer->SetXTitle("sftlayer");
  hsftlayer->Draw();
  
  TCanvas *csftseg = new TCanvas("csftseg","sft seg",800,0,400,400);
  csftseg->cd();
  hsftseg->SetXTitle("sft ch");
  hsftseg->Draw();
  
  TCanvas *csftp = new TCanvas("csftp","sftp",0,500,400,400);
  csftp->cd();
  hsftp->SetXTitle("beam momentum [GeV/c]");
  hsftp->Draw();
  
  TCanvas *csftposlx = new TCanvas("csftposlx","csftposlx",400,500,400,400);
  csftposlx->cd();
  hsftposlx->SetXTitle("local x position (mm)");
  hsftposlx->Draw();

  TCanvas *csftposly = new TCanvas("csftposly","csftposly",800,500,400,400);
  csftposly->cd();
  hsftposly->SetXTitle("local y position (mm)");
  hsftposly->Draw();

  TCanvas *csftposgx = new TCanvas("csftposgx","csftposgx",0,1000,400,400);
  csftposgx->cd();
  hsftposgx->SetXTitle("global x position (mm)");
  hsftposgx->Draw();

  TCanvas *csftposgy = new TCanvas("csftposgy","csftposgy",400,1000,400,400);
  csftposgy->cd();
  hsftposgy->SetXTitle("global y position (mm)");
  hsftposgy->Draw();

  TCanvas *csftposgz = new TCanvas("csftposgz","csftposgz",800,1000,400,400);
  csftposgz->cd();
  hsftposgz->SetXTitle("global z position (mm)");
  hsftposgz->Draw();


  
  //TString *dirname = filename;
  int len = strlen(filename);
  char ifile[len+1];
  strcpy(ifile,filename);
  strtok(ifile,"_");
  char *date = strtok(0,"v");  
  int version = atoi(strtok(0,"."));
  std::cout << "date " << date << std::endl;
  std::cout << "version " << version << std::endl;
  char dirname[256];
  sprintf(dirname,"%sv%d",date,version);
  gSystem->mkdir(dirname);
  gSystem->cd(dirname);
  csftnhits->SaveAs("sftnhits.pdf");
  csftlayer->SaveAs("sftlayer.pdf");
  csftseg->SaveAs("sftseg.pdf");
  csftp->SaveAs("sftp.pdf");
  csftposlx->SaveAs("sftposlx.pdf");
  csftposly->SaveAs("sftposly.pdf");
  csftposgx->SaveAs("sftposgx.pdf");
  csftposgy->SaveAs("sftposgy.pdf");
  csftposgz->SaveAs("sftposgz.pdf");

}
