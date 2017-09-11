//author: Hidemitsu Asano
//last update: June 24th, 2017
//analysis code for proto type BSFT read out by 2 NIM-EASIROCs
//
//This code is compiled by "make all" (only "make" is not sufficient)
//Main function of this code is to produce disk resident files.
//Histograms and TTrees to be looked should be defined in "InitializeHistograms()" 
//

#include "VEvent.hh"

#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

#include "ConfMan.hh"
#include "HistHelper.hh"
#include "DetectorID.hh"
#include "GeomMan.hh"
#include "TrParameters.hh"
#include "RawData.hh"
#include "TrHit.hh"
#include "SFTRawHit.hh"
#include "SFTCluster.hh"
#include "BSFTReco.hh"
#include "TrLocalTrack.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

//const int Verbosity = 100;
const int Verbosity = 0;
bool is_maketree = false;
bool is_tdccut = true;
double pede_mean[128];
double pede_sigma[128];
const char XUV[3][2] = {"X","U","V"};
const char XUVorder[12][2] = {"X","U","V","X","U","V","U","V","X","U","V","X"};

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class AnaBFT 
  : public VEvent
{

private:
  RawData *rawData;
  BSFTReco *BSFTRecoPtr;
  bool is_doTracking;
public:
  AnaBFT();
  ~AnaBFT();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  //bool InitializeHistograms(); 
  void InitializeEvent();
};

AnaBFT::AnaBFT()
  : VEvent(),
    rawData(0),
    BSFTRecoPtr(new BSFTReco()),
    is_doTracking(true)
{
  //InitializeHistograms(); 
}

AnaBFT::~AnaBFT()
{
  if (BSFTRecoPtr)   delete BSFTRecoPtr;
  if (rawData) delete rawData;
}

/*
struct FiberHit{
  int layer;
  int fiber;
  int easiroc;
  int type;
  double pos;
  int adchigh;
  bool otradchigh;
  int adclow;
  bool otradclow;
  int tdcleading;
  int tdctrailing;

  FiberHit(){
    layer = -1;
    fiber = -1;
    easiroc = -1;
    type = -1;
    pos = -99.;
    adchigh = -1;
    otradchigh = false;
    adclow = -1;
    otradclow = false;
    tdcleading = -1; 
    tdctrailing = -1;
  }
  void clear(){
    layer = -1;
    fiber = -1;
    easiroc = -1;
    type = -1;
    pos = -99.;
    adchigh = -1;
    otradchigh = false;
    adclow = -1;
    otradclow = false;
    tdcleading = -1; 
    tdctrailing = -1;
  }
};
static FiberHit sfiber;
*/


struct Event{
  unsigned int nevents;

  //BFT rawhit based info.
  std::vector<int> vTch;//0-128
  std::vector<int> vTlayer;//0-3
  std::vector<int> vTfiber;//0-12 ?
  std::vector<int> vTeasiroc;//0-4
  std::vector<int> vTtype;//0:X,U:1,V:2
  //std::vector<double> vTlpos;//local pos.
  std::vector<int> vTadchigh;//0-4095
  std::vector<bool> vTotradchigh;//0-1
  //std::vector<int> vTadclow;//0-4095
  //std::vector<bool> vTotradclow;//0-1
  std::vector<int> vTtdcLead;//0-4095
  std::vector<int> vTtdcTrail;//0-4095

  //BFT cluster based info.
  int    sftnhits;
  std::vector<int> sftlayer;
  std::vector<double> sftposx, sftposy;//always 0 if realistic detector is choosed
  int    sftnclus;
  std::vector<int> sftlayerc;
  std::vector<int> sftclssize;
  std::vector<double> sftclssizelx, sftclssizelz;//
  std::vector<double> sftposlxc, sftposgzc;//
  std::vector<double> sftposgxc,sftposgyc;//global position calculated from localx and tiltangle

  //Local tracking ( search beam position at z=0 in global coordinate system.
  int    ntr;
  std::vector<int>    layer;
  std::vector<double> chisqr;
  std::vector<double> x0, y0;
  std::vector<double> x0diff,y0diff;
  std::vector<double> u0, v0;
  std::vector<double> pos, res;
  std::vector<double> sftxangle,sftyangle;//track angle at layer 0
  std::vector<double> xanglediff,yanglediff;//track angle at layer 0
};
static Event event;


struct Cluster
{
  unsigned int nevents;

  int sftlayerc;
  int sftclssize;
  double sftclssizelx, sftclssizelz;//
  double sftposlxc, sftposgzc;//
  double sftposgxc,sftposgyc;//global position calculated from localx and tiltangle

};
static Cluster cluster;



bool AnaBFT::ProcessingBegin()
{
 //InitializeHistograms();   
 ConfMan *confMan=ConfMan::GetConfManager();
 if(!confMan){
   std::cout << "Can not find ConfManager !! " << std::endl;
   return false;
 }

 GeomMan *geomMan=confMan->GetGeomManager();
 if(!geomMan){
   std::cout << "Can not find GeomManager !! " << std::endl;
   return false;
 }
 
 for(int ich=0;ich<128;ich++){
   int type = geomMan->getXUV(ich);
   std::cout << "ch:" << ich << " layer:" << geomMan->getlayer(ich) 
        << " " << XUV[type] << " fiber:" 
        << geomMan->getfiber(ich) << std::endl;
 }


}

bool AnaBFT::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";
 

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cout << "Can not find ConfManager !! " << std::endl;
    return false;
  }
  
  GeomMan *geomMan=confMan->GetGeomManager();
  if(!geomMan){
    std::cout << "Can not find GeomManager !! " << std::endl;
    return false;
  }


  InitializeEvent();
  rawData = new RawData;
  if( !rawData->Decode(In) ){
    std::cerr << "event " << event.nevents  << " Decoding error! " << std::endl;
    std::cerr << "skipping this event.... " << std::endl;
    std::cerr << std::endl;
    
    return false;
  } 
  
  TTree *tree=NULL;
  if(is_maketree){
    tree = static_cast<TTree *>(gFile->Get("tree"));
  }//TTree *ctree = static_cast<TTree *>(gFile->Get("ctree"));
  
  const int totalch = 128;
  TH1I* adcHigh[totalch];
  TH1I* tdcLeading[totalch];
  TH1I* tdcTrailing[totalch];
  TH1I* tdcToT[totalch];
  TH2I* adcToTcorr[totalch];
 // int nbin = 4096;
  
  for(int ich = 0; ich < totalch; ++ich) {
    char hname[256];
    sprintf(hname,"ADC_HIGH_%d", ich);
    adcHigh[ich] = static_cast<TH1I*>(gFile->Get(hname)); 
    sprintf(hname,"TDC_LEADING_%d", ich);
    tdcLeading[ich] = static_cast<TH1I*>(gFile->Get(hname));  
    sprintf(hname,"TDC_TRAILING_%d", ich);
    tdcTrailing[ich] = static_cast<TH1I*>(gFile->Get(hname));  
    sprintf(hname,"TDC_TOT_%d",ich);
    tdcToT[ich] = static_cast<TH1I*>(gFile->Get(hname));
    sprintf(hname,"AdcToTcorr_%d",ich);
    adcToTcorr[ich] = static_cast<TH2I*>(gFile->Get(hname));
  }
  
  TH2I* rawhitcorr_ch[NumOfLayersSFT][3];
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    for(int ifirst=0;ifirst<3;ifirst++){
      char hname[256];
      sprintf(hname,"rawhitcorr%d_%d",ilr,ifirst);
      rawhitcorr_ch[ilr][ifirst] = static_cast<TH2I*>(gFile->Get(hname));
    }
  }

  TH1I* hitprofile_ch[NumOfLayersSFT];//layer , type(xuv)
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"Hitproflayer_ch%s%d",XUVorder[ilr],ilr);
    hitprofile_ch[ilr] = static_cast<TH1I*>(gFile->Get(hname));
  }
  
  TH1I* hitmulti[NumOfLayersSFT];//
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"HitMulti%s%d",XUVorder[ilr],ilr);
    hitmulti[ilr] = static_cast<TH1I*>(gFile->Get(hname));
  }

  TH2F* hX0corr[NumOfLayersSFT];
  for(int ilr=0;ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"hX0corr%d",ilr);
    hX0corr[ilr] = static_cast<TH2F*>(gFile->Get(hname));
  }

  TH1I* hclsadc[NumOfLayersSFT];
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"hclsadc%d",ilr);
    hclsadc[ilr] = static_cast<TH1I*>(gFile->Get(hname));
  }

  //SFT rawhit info.
  {
    const SFTRawHitContainer &cont = rawData->GetSFTRawHitContainer();
    int nh=cont.size();
    //std::cout << "N raw hit " << nh << std::endl;
    event.sftnhits = nh;
    if(nh!=128){
      std::cout << __FILE__ << " l." << __LINE__ << "size of container " << nh << std::endl; 
    }
    for( int ich=0; ich<nh; ++ich ){
      SFTRawHit *hit=cont[ich];
      //event.sftlayer.push_back(hit->LayerId());
      event.vTlayer.push_back(hit->LayerId());
      event.vTch.push_back(hit->ChId());
      int adcval = hit->GetAdcHigh();
      event.vTadchigh.push_back(adcval);
      adcHigh[ich]->Fill(adcval);
      //event.vTotradchigh.push_back(hit->GetotrAdcHigh());
      //event.vTadclow.push_back(hit->GetAdcLow());
      //event.vTotradclow.push_back(hit->GetotrAdcLow());
      int tdclval= hit->GetTdcLeading();
      event.vTtdcLead.push_back(tdclval);
      tdcLeading[ich]->Fill(tdclval);
      int tdctval= hit->GetTdcTrailing();
      event.vTtdcTrail.push_back(tdctval);
      tdcTrailing[ich]->Fill(tdctval);
      int totval = tdclval - tdctval;
      tdcToT[ich]->Fill(totval);
      adcToTcorr[ich]->Fill(adcval,totval);
      //int layer = hit->LayerId();
      //int fiber = hit->FiberId();
      //int type = hit->GetType();
    }
  }

  
  BSFTRecoPtr->SetSigmaThreshold (10.);
  BSFTRecoPtr->UseTOTcut(false);
  BSFTRecoPtr->MakeSFTRawHits( rawData ); // name of this func. is confusing in real data analysis. change it
  //SFT TrHit (hits above threshold)
  //TODO change name of class "TrHit", these hits are not necessarily associated track.
   
  const int maxmul=10;
  int multi[NumOfLayersSFT]={-1};
  int fibercont[NumOfLayersSFT][maxmul];
  for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
    for(int imul =0;imul<maxmul;imul++){
      fibercont[ilr][imul]=-1;
    }
  }

  for( int ilr=0; ilr<NumOfLayersSFT; ++ilr ){
    const TrHitContainer trhitcont = BSFTRecoPtr->GetSFTTrHitContainer(ilr);
    int nhit = trhitcont.size();
    hitmulti[ilr]->Fill(nhit);
    multi[ilr] = nhit;
    for(int ihit=0;ihit<nhit ; ihit++){
      TrHit *hit = trhitcont[ihit];
      //std::cout << std::endl;
      //hit->Print();
      int layer = hit->GetLayer();
      if(layer != ilr){
        std::cout << "layer number does not match !! " << std::endl;
        return false;
      }
      int fiber = hit->GetFiber();
      hitprofile_ch[ilr]->Fill(fiber);
      fibercont[ilr][ihit]=fiber;      
      if(fiber<0) hit->Print();
    }
  }

  for( int ilr=0; ilr<NumOfLayersSFT; ++ilr ){
    for( int ifirst=0; ifirst<3; ++ifirst){
      for( int imul0=0; imul0<multi[ifirst]; ++imul0){
        for( int imul=0; imul<multi[ilr]; ++imul){
          rawhitcorr_ch[ilr][ifirst]->Fill(fibercont[ifirst][imul0],fibercont[ilr][imul]);
        }
      }
    }
  }



  BSFTRecoPtr->SortTrHits();//sort TrHits by segment ID
  //TODO check the TrHits are sorted or not, here.
  BSFTRecoPtr->SFTClustering();//clustering, 
  
  TH1F* hitprofile[NumOfLayersSFT];//
  TH1I* hitclsmulti[NumOfLayersSFT]={NULL};
  TH1I* hitsize[NumOfLayersSFT]={NULL};
  
  for(int ilr =0; ilr<NumOfLayersSFT; ilr++){
    char hname[256];
    sprintf(hname,"Hitproflayer%s%d",XUVorder[ilr],ilr);
    hitprofile[ilr] = static_cast<TH1F*>(gFile->Get(hname));

    sprintf(hname,"HitClsMulti%s%d",XUVorder[ilr],ilr);
    hitclsmulti[ilr] = static_cast<TH1I*>(gFile->Get(hname));
    
    sprintf(hname,"HitSize%s%d",XUVorder[ilr],ilr);
    hitsize[ilr] = static_cast<TH1I*>(gFile->Get(hname));
  }

  int nclsX0 = 0;
  double globalx0=0;
  for( int ilr=0; ilr<NumOfLayersSFT; ++ilr ){
    const SFTClusterContainer cluscont = BSFTRecoPtr->GetSFTClusterContainer(ilr);
    int nclus = cluscont.size();
    if(ilr==0) nclsX0 = nclus;
    event.sftnclus = nclus;
    if(nclus!=0) hitclsmulti[ilr]->Fill(nclus);
    if(Verbosity>2){
      std::cout << __FILE__ << "  " << __LINE__ << " nclus " << nclus << std::endl;
    }
    for( int iclus=0; iclus<nclus;iclus++){
      SFTCluster *clus = cluscont[iclus];
      int layer = clus->GetLayer();
      if(layer != ilr){
        std::cout << "invalid layer ID !! something is wrong" << std::endl;
        return -1;
      }
      unsigned int clsID = clus->GetClusterID();
      int clustersize = clus->GetClusterSize();
      hitsize[ilr]->Fill(clustersize);
      double clusterlxsize = clus->GetClusterLxSize();
      double clusterlzsize = clus->GetClusterLzSize();
      double lx = clus->GetLocalX();
      hitprofile[ilr]->Fill(lx);
      double gz = clus->GetGlobalZ();
      double angle = clus->GetTiltAngle()*Deg2Rad;//degree to radian
      double globalx = lx * cos(angle);
      double globaly = lx * sin(angle);
      if( (ilr==0) && (nclsX0>0)) globalx0 = globalx;
     
     if( nclsX0>0 ){
        hX0corr[ilr]->Fill(globalx0,globalx);
        //std::cout << "angle " << angle << std::endl;
        //std::cout << "cos(angle) " <<  cos(angle) << std::endl;
        //std::cout << "lx         " << lx  << std::endl;
        //std::cout << "globalx    " << globalx << std::endl;
        //std::cout << "X0         " << globalx0 << std::endl;
      }
      int adcsum = clus->GetAdcSum();
      hclsadc[ilr]->Fill(adcsum);

      if(Verbosity>3){
        std::cout << "layer " << ilr << std::endl;
        std::cout << "clsID " << clsID << std::endl;
        std::cout << "lx " << lx << std::endl;
        std::cout << "global z " << gz << std::endl;
      }
      cluster.nevents = event.nevents;
      //double adc = clus->GetAdcSum();
      event.sftlayerc.push_back(ilr);
      event.sftclssize.push_back(clustersize);
      event.sftclssizelx.push_back(clusterlxsize);
      event.sftclssizelz.push_back(clusterlzsize);
      event.sftposlxc.push_back(lx);
      event.sftposgzc.push_back(gz);
      event.sftposgxc.push_back(globalx);
      event.sftposgyc.push_back(globaly);

      cluster.sftlayerc = ilr;
      cluster.sftclssize = clustersize;
      cluster.sftclssizelx = clusterlxsize;
      cluster.sftclssizelz = clusterlzsize;
      cluster.sftposlxc = lx;
      cluster.sftposgzc = gz;
      cluster.sftposgxc = globalx;
      cluster.sftposgyc = globaly;
      //ctree->Fill();

      if(Verbosity>3){
        std::cout << __FILE__ << "  " << __LINE__ << " event  " << event.nevents << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " cluster  " << iclus << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " layer  " << ilr << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " lx  " << lx << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " global z  " << gz << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " global x  " << globalx << std::endl;
        std::cout << __FILE__ << "  " << __LINE__ << " global y  " << globaly << std::endl;
      }
    }//for iclus
  }//for ilr

    
  TH2F* hxy0 = static_cast<TH2F*>(gFile->Get("hxy0"));
  TH1I* hnassociate = static_cast<TH1I*>(gFile->Get("hnassociate"));
  TH1F* hchi2 = static_cast<TH1F*>(gFile->Get("hchi2"));
  TH1F* hresi[NumOfLayersSFT];
  for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
    char hname[256];
    sprintf(hname,"hresi%s%d",XUVorder[ilr],ilr);
    hresi[ilr] = static_cast<TH1F*>(gFile->Get(hname));
  }
  


  if(is_doTracking){  
    BSFTRecoPtr->TrackSearchSFTT();//making index
    //This loop is for evaluating tracking performance while event matching between
    //primary info and reconstructed track info.
    int nt=BSFTRecoPtr->GetNtracksSFTT();
    event.ntr=nt;
    if(Verbosity>10){
      std::cout << __FILE__ << "  " << __LINE__ << " ntrack  " << nt << std::endl;
    }
    for( int it=0; it<nt; ++it ){
      TrLocalTrack *tp=BSFTRecoPtr->GetTrackSFTT(it);
      int nh=tp->GetNHit();
      hnassociate->Fill(nh);
      //get chi2/ndf
      double chisqr=tp->GetChiSquare();
      hchi2->Fill(chisqr);
      //double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();//fitting results of u0 and v0 position
                                            //TODO : what is the definition of u0 and v0 ?
      double xtgt=tp->GetX( 0.0 ), ytgt=tp->GetY( 0.0 );
      double utgt=u0, vtgt=v0;
      hxy0->Fill(xtgt,ytgt);

      event.chisqr.push_back(chisqr);
      event.x0.push_back(xtgt);
      event.y0.push_back(ytgt);
      event.u0.push_back(utgt);
      event.v0.push_back(vtgt); 


      
      for( int ih=0; ih<nh; ++ih ){
        SFTCluster *hit=tp->GetHit(ih);
        int layerId=hit->GetLayer(); 
        event.layer.push_back(layerId);  
        double pos=hit->GetLocalX(), res=hit->GetResidual();
        hresi[layerId]->Fill(res);
        event.pos.push_back(pos);
        event.res.push_back(res);
        
        //check angle resolution at layer 0
        if(layerId==0){
          double projxL0 = hit->GetProjectedPosX();
          double projyL0 = hit->GetProjectedPosY();
          double gzL0 = hit->GetGlobalZ();
          double xangle = atan2(projxL0-xtgt,(double)gzL0);
          double yangle = atan2(projyL0-ytgt,gzL0);
      //    double xangle = atan2(projxL0,gzL0);
      //    double yangle = atan2(projyL0,gzL0);
          
          if(Verbosity>0){
            std::cout << "EVENT " << event.nevents  << std::endl;
            std::cout << "projxL0 " << projxL0  << std::endl;
            std::cout << "xtgt " << xtgt  << std::endl;
            std::cout << "diff " << projxL0-xtgt  << std::endl;
            std::cout << "atan2 " << xangle  << std::endl;
          }
          event.sftxangle.push_back(xangle);
          event.sftyangle.push_back(yangle);
        }
      }
    }
  }//end of tracking

  event.nevents++;
  if(is_maketree) tree->Fill();
  if(Verbosity>0){
    std::cout << __FILE__ << "  " << __LINE__ << " fill " << std::endl;
  }
  return true;
}

void AnaBFT::InitializeEvent( void )
{
  //SFT rawhit
  event.vTch.clear();
  event.vTeasiroc.clear();
  event.vTtype.clear();
  //event.vTlpos.clear();
  event.vTadchigh.clear();
  event.vTotradchigh.clear();
  //event.vTadclow.clear();
  //event.vTotradclow.clear();
  event.vTtdcLead.clear();
  event.vTtdcTrail.clear();


  //SFT cluser
  event.sftnhits = 0;
  event.sftlayer.clear();
  event.sftposx.clear();
  event.sftposy.clear();
  event.sftnclus = 0;
  event.sftlayerc.clear();
  event.sftclssize.clear();
  event.sftclssizelx.clear();
  event.sftclssizelz.clear();
  event.sftposlxc.clear();
  event.sftposgzc.clear();
  event.sftposgxc.clear();
  event.sftposgyc.clear();
  event.sftxangle.clear();
  event.sftyangle.clear();

  //SFT Local Tracking
  event.ntr = -1;
  event.layer.clear();
  event.chisqr.clear();
  event.x0.clear();
  event.y0.clear();
  event.x0diff.clear();
  event.y0diff.clear();
  event.u0.clear();
  event.v0.clear();
  event.pos.clear();
  event.res.clear();
  event.xanglediff.clear();
  event.yanglediff.clear();
  
  if(Verbosity>0){
    std::cout << __FILE__ << "  " << __LINE__ << " Init. end " << std::endl;
  }
}

bool AnaBFT::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new AnaBFT;
}

bool ConfMan::InitializeHistograms()
{ 
  
  //SFT cluster based tree
  if(is_maketree){
    HBTree("ctree","cluster based tree");
    TTree *ctree = static_cast<TTree *>(gFile->Get("ctree"));
    ctree->Branch("event", &cluster.nevents);
    ctree->Branch("sftlayerc", &cluster.sftlayerc);
    ctree->Branch("sftclssize", &cluster.sftclssize);//cluster size = number of hit segment
    ctree->Branch("sftclssizelx", &cluster.sftclssizelx);//cluster size in local x coordinate
    ctree->Branch("sftclssizelz", &cluster.sftclssizelz);//cluster size in local z coordiante
    ctree->Branch("sftposlxc", &cluster.sftposlxc);//local hit position of cluster
    ctree->Branch("sftposgzc", &cluster.sftposgzc);//global hit position of cluster
    ctree->Branch("sftposgxc", &cluster.sftposgxc);//global hit position of cluster
    ctree->Branch("sftposgyc", &cluster.sftposgyc);//global hit position of cluster

    HBTree("tree","event wise tree");
    TTree *tree = static_cast<TTree *>(gFile->Get("tree"));
    //tree = new TTree("tree","event wise tree");
    tree->Branch("event", &event.nevents);
    //SFT rawhit
    tree->Branch("ch",&event.vTch);
    tree->Branch("layer",&event.vTlayer);
    tree->Branch("fiber",&event.vTfiber);
    tree->Branch("easiroc",&event.vTeasiroc);
    tree->Branch("type",&event.vTtype);
    //tree->Branch("lpos",&event.vTlpos);
    tree->Branch("adchigh",&event.vTadchigh);
    tree->Branch("otradchigh",&event.vTotradchigh);
    //tree->Branch("adclow",&event.vTadclow);
    //tree->Branch("otradclow",&event.vTotradclow);
    tree->Branch("tdclead",&event.vTtdcLead);
    tree->Branch("tdctrail",&event.vTtdcTrail);
    
    //SFT cluster
    tree->Branch("sftnhits", &event.sftnhits);
    tree->Branch("sftlayer", &event.sftlayer);
    tree->Branch("sftposx", &event.sftposx);//hit position of raw hits
    tree->Branch("sftposy", &event.sftposy);//hit position of raw hits
    tree->Branch("sftlayerc", &event.sftlayerc);
    tree->Branch("sftclssize", &event.sftclssize);//cluster size = number of hit segment
    tree->Branch("sftclssizelx", &event.sftclssizelx);//cluster size in local x coordinate
    tree->Branch("sftclssizelz", &event.sftclssizelz);//cluster size in local z coordiante
    tree->Branch("sftposlxc", &event.sftposlxc);//local hit position of cluster
    tree->Branch("sftposgzc", &event.sftposgzc);//global hit position of cluster
    tree->Branch("sftposgxc", &event.sftposgxc);//global hit position of cluster
    tree->Branch("sftposgyc", &event.sftposgyc);//global hit position of cluster

    //Local Tracking
    tree->Branch("ntr", &event.ntr);
    tree->Branch("layer", &event.layer);
    tree->Branch("chisqr", &event.chisqr);
    tree->Branch("x0", &event.x0);
    tree->Branch("y0", &event.y0);
    tree->Branch("x0diff", &event.x0diff);
    tree->Branch("y0diff", &event.y0diff);
    tree->Branch("u0", &event.u0);
    tree->Branch("v0", &event.v0);
    tree->Branch("pos", &event.pos);
    tree->Branch("res", &event.res);
    tree->Branch("sftxangle", &event.sftxangle);//x angle
    tree->Branch("sftyangle", &event.sftyangle);//y angle
    tree->Branch("xanglediff",&event.xanglediff);
    tree->Branch("yanglediff",&event.yanglediff);
  }
  
  //raw infomation
  const int totalch = 128;
  TH1I* adcHigh[totalch];
  TH1I* tdcLeading[totalch];
  TH1I* tdcTrailing[totalch];
  TH1I* tdcToT[totalch];
  TH2I* adcToTcorr[totalch];
  int nbin = 4096;
  for(int ich = 0; ich < totalch; ++ich) {
    adcHigh[ich] = new TH1I(Form("ADC_HIGH_%d", ich),
        Form("ADC high gain %d", ich),
        nbin, 0, 4096);
    adcHigh[ich]->SetXTitle("ch.");
    tdcLeading[ich] = new TH1I(Form("TDC_LEADING_%d", ich),
        Form("TDC leading %d", ich),
        nbin, 0, 4096);
    tdcLeading[ich]->SetXTitle("ch.");
    tdcTrailing[ich] = new TH1I(Form("TDC_TRAILING_%d", ich),
        Form("TDC trailing %d", ich),
        nbin, 0, 4096);
    tdcTrailing[ich]->SetXTitle("ch.");
    tdcToT[ich] = new TH1I(Form("TDC_TOT_%d",ich),
      Form("Time Over Threshold %d",ich),
      512,0,512);
    tdcToT[ich]->SetXTitle("ch.");
    adcToTcorr[ich] = new TH2I(Form("AdcToTcorr_%d",ich),
      Form("ADC vs TOT correlation %d", ich),
        256,0,4096,256,0,256);
    adcToTcorr[ich]->SetXTitle("ADC [ch]");
    adcToTcorr[ich]->SetYTitle("TOT [ch]");
  }

  for(int ich = 0; ich < totalch; ich++){
    pede_mean[ich] = pedestalmean_[ich];
    pede_sigma[ich] = pedestalsigma_[ich];
  }


  TH1I* hitprofile_ch[NumOfLayersSFT];//
  for(int ilr=0 ;ilr<NumOfLayersSFT; ilr++){
    hitprofile_ch[ilr] = new TH1I(Form("Hitproflayer_ch%s%d",XUVorder[ilr],ilr),Form("Hit ch. profile %s layer %d",XUVorder[ilr],ilr), 12,0,12); 
    hitprofile_ch[ilr]->SetXTitle("fiber ch.#");
  }

  //raw hit mulitiplicity
  TH1I* hitmulti[NumOfLayersSFT];
  for(int ilr =0 ; ilr<NumOfLayersSFT; ilr++){
    hitmulti[ilr] = new TH1I(Form("HitMulti%s%d",XUVorder[ilr],ilr),Form("Hit Multiplicity %s layer %d",XUVorder[ilr],ilr), 12,0,12); 
    hitmulti[ilr]->SetXTitle("hit multiplicity");
  }
  
  TH2I* rawhitcorr_ch[NumOfLayersSFT][3];
  for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
    for(int ifirst = 0 ; ifirst<3 ; ifirst++){
      rawhitcorr_ch[ilr][ifirst] = new TH2I(Form("rawhitcorr%d_%d",ilr,ifirst),
      Form("raw Hit correlation %s layer %d - %s layer %d",XUVorder[ifirst],ifirst,XUVorder[ilr],ilr),13,0,13,13,0,13);
      rawhitcorr_ch[ilr][ifirst]->SetXTitle("fiber #");
      rawhitcorr_ch[ilr][ifirst]->SetYTitle("fiber #");
    }
  }

  //cluster information
  TH1F* hitprofile[NumOfLayersSFT];//
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    hitprofile[ilr] = new TH1F(Form("Hitproflayer%s%d",XUVorder[ilr],ilr),Form("Hit profile %s layer %d after clustering",XUVorder[ilr],ilr), 40,-10,10); 
    hitprofile[ilr]->SetXTitle("local pos. [mm]");
  }
  
  TH1I* hitclsmulti[NumOfLayersSFT];//
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    hitclsmulti[ilr] = new TH1I(Form("HitClsMulti%s%d",XUVorder[ilr],ilr),Form("Hit Multiplicity %s layer %d after clustering",XUVorder[ilr],ilr), 20,0,20); 
    hitclsmulti[ilr]->SetXTitle("hit multiplicity");
  }
  
  TH1I* hitsize[NumOfLayersSFT];//cluster size
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    hitsize[ilr] = new TH1I(Form("HitSize%s%d",XUVorder[ilr],ilr),Form("cluster size %s layer %d",XUVorder[ilr],ilr), 12,0,12); 
    hitsize[ilr]->SetXTitle("cluster size");
  }

  //cluste position correlation btw X0 - other layers
  TH2F* hX0corr[NumOfLayersSFT];
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){  
    hX0corr[ilr] = new TH2F(Form("hX0corr%d",ilr),Form("Hit correlation X0 - %s layer %d after clustering",XUVorder[ilr],ilr), 40,-4,4,40,-4,4); 
    hX0corr[ilr]->SetXTitle("X0 pos. [mm]");
    hX0corr[ilr]->SetYTitle("pos. [mm]");
  }

  TH1I* hclsadc[NumOfLayersSFT];
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    hclsadc[ilr] = new TH1I(Form("hclsadc%d",ilr),Form("ADC sum distibution layer %d %s", ilr,XUVorder[ilr]),1000,0,10000);  
    hclsadc[ilr]->SetXTitle("ADC sum [ch]");
  }

  
  TH1I* hnassociate = new TH1I("hnassociate","number of associated hits",15,0,15);
  hnassociate->SetXTitle("nhit");
  
  TH2F* hxy0 = new TH2F("hxy0","X0:Y0 pos.", 300,-30,30,300,-30,30);
  hxy0->GetXaxis()->SetTitle("x0 [mm]");
  hxy0->GetYaxis()->SetTitle("y0 [mm]");

  TH1F* hchi2 = new TH1F("hchi2","chi2/ndf", 1000,0,200);
  hchi2->SetXTitle("chi2/ndf");

  TH1F* hresi[NumOfLayersSFT];
  for(int ilr=0; ilr<NumOfLayersSFT; ilr++){
    hresi[ilr] = new TH1F(Form("hresi%s%d",XUVorder[ilr],ilr),Form("residual %s layer %d",XUVorder[ilr],ilr),100,-10,10);
    hresi[ilr]->SetXTitle("residual [mm]");
  }


  std::cout << "Init "<< std::endl;
  return true;

}
