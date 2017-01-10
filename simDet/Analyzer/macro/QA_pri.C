//QA_pri
//quick check for primary information and generate pdf files
//


void QA_pri(const char* filename)
{
  gStyle->SetOptStat("emr");
  //TFile::Open(filename);
  TFile *tf = new TFile(filename, "READ");
  TTree *tree = (TTree*)(tf->Get("tree"));
  
  double bpriposx=0;
  double bpriposy=0;
  double bpriposz=0;
  double bpbeam=0;
  double bubeam=0;
  double bvbeam=0;
  double bxangle=0;
  double byangle=0;
  
  tree->SetBranchAddress("priposx",&bpriposx);
  tree->SetBranchAddress("priposy",&bpriposy);
  tree->SetBranchAddress("priposz",&bpriposz);
  tree->SetBranchAddress("pbeam",&bpbeam);
  tree->SetBranchAddress("ubeam",&bubeam);
  tree->SetBranchAddress("vbeam",&bvbeam);
  tree->SetBranchAddress("xangle",&bxangle);
  tree->SetBranchAddress("yangle",&byangle);


  TH1D *hpriposx = new TH1D("hpriposx","primary x pos.",600,-60,60);
  TH1D *hpriposy = new TH1D("hpriposy","primary y pos.",600,-60,60);
  TH1D *hpriposz = new TH1D("hpriposz","primary z pos.",600,-60,60);
  TH1D *hpbeam   = new TH1D("hpbeam","primary beam momentum ",200,19.5,20.5);
  TH1D *hubeam   = new TH1D("hubeam","differential primary beam momentum dP_{x}/dP_{z}",200,-0.02,0.02);
  TH1D *hvbeam   = new TH1D("hvbeam","differential primary beam momentum dP_{y}/dP_{z}",200,-0.02,0.02);
  TH1D *hxangle  = new TH1D("hxangle","primary beam angle in x",200,-0.02,0.02);
  TH1D *hyangle  = new TH1D("hyangle","primary beam angle in y",200,-0.02,0.02);

  Long64_t nentries = tree->GetEntries();
  std::cout << "# of entries " << nentries << std::endl;
  for(Long64_t i=0;i < nentries;i++){
    tree->GetEntry(i);
    hpriposx->Fill(bpriposx); 
    hpriposy->Fill(bpriposy); 
    hpriposz->Fill(bpriposz); 
    hpbeam->Fill(bpbeam);
    hubeam->Fill(bubeam);
    hvbeam->Fill(bvbeam);
    hxangle->Fill(bxangle);
    hyangle->Fill(byangle);
  }
  


  TCanvas *cpriposx = new TCanvas("cpriposx","primary pos x",0,0,400,400);
  cpriposx->cd();
  hpriposx->SetXTitle("position x [mm]");
  hpriposx->Draw();

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


}
