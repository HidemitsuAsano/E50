void DrawClusterInfo(const char *filename="20170615T144153_000015_0.root"){


  TFile *fin = new TFile(filename,"READ");
  if(!fin) return;
  gStyle->SetOptStat("e");
  gStyle->SetStatY(0.9);
  //cluster size
  TCanvas *cmulti = new TCanvas("cmulti","cmulti",1200,800);
  TCanvas *csize = new TCanvas("csize","cluster size",1200,800);
  cmulti->Divide(4,3);
  csize->Divide(4,3);
  
  const int NumOfLayersSFT = 12;
  const char XUVorder[12][2] = {"X","U","V","X","U","V","U","V","X","U","V","X"};
  TH1I* hitclsmulti[NumOfLayersSFT]={NULL};
  TH1I* hitsize[NumOfLayersSFT]={NULL};
  for(int ilr =0; ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"HitClsMulti%s%d",XUVorder[ilr],ilr);
    hitclsmulti[ilr] = (TH1I*)(fin->Get(hname));
    cmulti->cd(ilr+1);
    hitclsmulti[ilr]->Draw();

    sprintf(hname,"HitSize%s%d",XUVorder[ilr],ilr);
    hitsize[ilr] = (TH1I*)(fin->Get(hname));
    csize->cd(ilr+1);
    hitsize[ilr]->Draw();
  }

  TCanvas *chitprof = new TCanvas("chitprof","chitprof",1200,800);
  chitprof->Divide(4,3);
  TH1F* hitprofile[NumOfLayersSFT];//
  for(int ilr =0; ilr<NumOfLayersSFT; ilr++){
    chitprof->cd(ilr+1);
    char hname[256];
    sprintf(hname,"Hitproflayer%s%d",XUVorder[ilr],ilr);
    hitprofile[ilr] = (TH1F)(gFile->Get(hname));
    hitprofile[ilr]->Draw();
  }
}

