//QA_gensft
//quick check for BSFT in GEANT4 info. and generate pdf files
#include <iostream>
#include <math.h>
#include <map>

#include <TRint.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>

void QA_gensft(const char* filename)
{
  gStyle->SetOptStat("emr");
  TFile *tf = new TFile(filename, "READ");
  TTree *tree = (TTree*)(tf->Get("tree"));
  
  //number of hits/event
  int  bsftnhits;
  //hit layer /hit segment (=fiber/MPPC)
  std::vector<int> *bsftlayer,*bsftseg;
  //PID on SFT, track No. on SFT
  std::vector<int> *bsftpid,*bsftid;
  //E deposit on SFT and time ?
  std::vector<double> *bsftedep, *bsfttime;
  //absolute value of Momentum  
  std::vector<double> *bsftp;
  //theta and phi of momentum vector of beam
  std::vector<double> *bstdtheta,*bsftphi;
  //local x and y position of sft 
  std::vector<double> *bsftposlx, *bsftposly;
  //global x, y and z position of sft
  std::vector<double> *bsftposgx, *bsftposgy, *bsftposgz;
  //pathlength , beta of particle
  std::vector<double> *bsftpath, *bsftbeta;
  //mass of particle
  std::vector<double> *bsftmass;
  
 
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
  std::cout << "# of entries " << nentries << std::endl;
  for(Long64_t i=0;i < nentries;i++){
    tree->GetEntry(i);
    hsftnhits->Fill(bsftnhits);
    if(bsftnhits != (int)bsftlayer->size()){
      std::cout << "something is wrong" << std::endl;
    }
    if(bsftnhits != (int)bsftseg->size()){
      std::cout << "something is wrong" << std::endl;
    }

    for(unsigned int ihit=0; ihit< bsftlayer->size(); ihit++){
      int sfthitlayer = bsftlayer->at(ihit);
      hsftlayer->Fill(sfthitlayer);
      int sfthitsegment = bsftseg->at(ihit);
      hsftseg->Fill(sfthitsegment);
    }
     
    hsftp->Fill(bsftp);
    hsftposlx->Fill(bsftposlx);
    hsftposly->Fill(bsftposly);
    hsftposgx->Fill(bsftposgx);
    hsftposgy->Fill(bsftposgy);
    hsftposgz->Fill(bsftposgz);


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
  
  TCanvas *csftposlx = new TCanvas("csftposlx",400,500,400,400);
  csftposlx->cd();
  hsftposlx->SetXTitle("local x position (mm)");

  TCanvas *csftposly = new TCanvas("csftposly",800,500,400,400);
  csftposly->cd();
  hsftposly->SetXTitle("local y position (mm)");

  TCanvas *csftposgx = new TCanvas("csftposgx",0,1000,400,400);
  csftposlx->cd();
  hsftposlx->SetXTitle("global x position (mm)");

  TCanvas *csftposgy = new TCanvas("csftposgy",400,1000,400,400);
  csftposgy->cd();
  hsftposgy->SetXTitle("global y position (mm)");

  TCanvas *csftposgz = new TCanvas("csftposgz",800,1000,400,400);
  csftposgz->cd();
  hsftposgz->SetXTitle("global z position (mm)");



  /*
  TCanvas *cpriposy = new TCanvas("cpriposy","primary pos y",400,0,400,400);
  cpriposy->cd();
  hpriposy->SetXTitle("position y [mm]");
  hpriposy->Draw();

  TCanvas *cpriposz = new TCanvas("cpriposz","primary pos z",800,0,400,400);
  cpriposz->cd();
  hpriposz->SetXTitle("position z [mm]");
  hpriposz->Draw();
  
  TCanvas *cpbeam = new TCanvas("cpbeam","primary beam momentum",0,500,400,400);
  cpbeam->cd();
  hpbeam->SetXTitle("beam momentum [GeV/c]");
  hpbeam->Draw();

  TCanvas *cubeam = new TCanvas("cubeam","cubeam",400,500,400,400);
  cubeam->cd();
  hubeam->SetXTitle("dP_{x}/dP_{z}");
  hubeam->Draw();

  TCanvas *cvbeam = new TCanvas("cvbeam","cvbeam",800,500,400,400);
  cvbeam->cd();
  hvbeam->SetXTitle("dP_{y}/dP_{z}");
  hvbeam->Draw();
  
  TCanvas *cxangle = new TCanvas("cxangle","cxangle",1200,500,400,400);
  cxangle->cd();
  hxangle->SetXTitle("primary beam angle in x [radian]");
  hxangle->Draw();

  TCanvas *cyangle = new TCanvas("cyangle","cyangle",1600,500,400,400);
  cyangle->cd();
  hyangle->SetXTitle("primary beam angle in y [radian]");
  hyangle->Draw();
  
  //TString *dirname = filename;
  char *ifile[strlen(filename)+1];
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
  cpriposx->SaveAs("priposx.pdf");
  cpriposy->SaveAs("priposy.pdf");
  cpriposz->SaveAs("priposz.pdf");
  cpbeam->SaveAs("pbeam.pdf");
  cubeam->SaveAs("ubeam.pdf");
  cvbeam->SaveAs("vbeam.pdf");
  cxangle->SaveAs("xangle.pdf");
  cyangle->SaveAs("yangle.pdf");
  */

}
