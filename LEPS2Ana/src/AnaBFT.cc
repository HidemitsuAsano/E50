//author: Hidemitsu Asano
//data: June 24th, 2017
//analysis code for proto type BFT readout by 2 NIM-EASIROCs

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
#include "RawData.hh"
#include "TrHit.hh"
#include "SFTCluster.hh"
#include "TrAnalyzer.hh"
#include "SpecLib.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

//const int Verbosity = 100;
const int Verbosity = 0;
bool is_maketree = false;
double pede_mean[128];
double pede_sigma[128];

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
  TrAnalyzer *TrAna;
  //  HodoAnalyzer *HodoAna;
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
    TrAna(new TrAnalyzer()),
    is_doTracking(true)
    //HodoAna(new HodoAnalyzer())
{
  //InitializeHistograms(); 
}

AnaBFT::~AnaBFT()
{
  //  if (HodoAna)   delete HodoAna;
  if (TrAna)   delete TrAna;
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
 std::cout << __FUNCTION__ << std::endl;
 return true;
}

bool AnaBFT::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  
  //ConfMan *confMan=ConfMan::GetConfManager();
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
     tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  }//TTree *ctree = dynamic_cast<TTree *>(gFile->Get("ctree"));
  
  const int totalch = 128;
  TH1I* adcHigh[totalch];
  TH1I* tdcLeading[totalch];
  TH1I* tdcTrailing[totalch];
  
 // int nbin = 4096;
  
  for(int i = 0; i < totalch; ++i) {
    char hname[256];
    sprintf(hname,"ADC_HIGH_%d", i);
    adcHigh[i] = dynamic_cast<TH1I*>(gFile->Get(hname)); 
    sprintf(hname,"TDC_LEADING_%d", i);
    tdcLeading[i] = dynamic_cast<TH1I*>(gFile->Get(hname));  
    sprintf(hname,"TDC_TRAILING_%d", i);
    tdcTrailing[i] = dynamic_cast<TH1I*>(gFile->Get(hname));  
  }

  //SFT
  {
    const SFTRawHitContainer &cont =rawData->GetSFTRawHitContainer();
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
      int val = hit->GetAdcHigh();
      event.vTadchigh.push_back(val);
      adcHigh[ich]->Fill(val);
      //event.vTotradchigh.push_back(hit->GetotrAdcHigh());
      //event.vTadclow.push_back(hit->GetAdcLow());
      //event.vTotradclow.push_back(hit->GetotrAdcLow());
      val= hit->GetTdcLeading();
      event.vTtdcLead.push_back(val);
      tdcLeading[ich]->Fill(val);
      val= hit->GetTdcTrailing();
      event.vTtdcTrail.push_back(val);
      tdcTrailing[ich]->Fill(val);
    }
  }
  
  
  //////////////////////////Tracking
  TrAna->SetSigmaThreshold (10.);
  TrAna->DecodeSFTRawHits( rawData ); // name of this func. is confusing in real data analysis. change it
  TrAna->SortTrHits();//sort TrHits by segment ID
  //TODO check the TrHits are sorted or not, here.
  TrAna->SFTClustering();//clustering, 

  for( int layer=0; layer<NumOfLayersSFT; ++layer ){
    const SFTClusterContainer cluscont = TrAna->GetSFTClusterContainer(layer);
    int nclus = cluscont.size();
    event.sftnclus = nclus;
    if(Verbosity>2){
      std::cout << __FILE__ << "  " << __LINE__ << " nclus " << nclus << std::endl;
    }
    for( int iclus=0; iclus<nclus;iclus++){
      SFTCluster *trhitclus = cluscont[iclus];
      unsigned int clsID = trhitclus->GetClusterID();
      int clustersize = trhitclus->GetClusterSize();
      double clusterlxsize = trhitclus->GetClusterLxSize();
      double clusterlzsize = trhitclus->GetClusterLzSize();
      double lx = trhitclus->GetLocalX();
      double gz = trhitclus->GetGlobalZ();
      double angle = trhitclus->GetTiltAngle();
      double globalx = lx * cos(angle);
      double globaly = lx * sin(angle);
    }
    if(Verbosity>3){
      std::cout << "layer " << layer << std::endl;
      std::cout << "clsID " << clsID << std::endl;
      std::cout << "lx " << lx << std::endl;
      std::cout << "global z " << gz << std::endl;
    }
    cluster.nevents = event.nevents;
    //double adc = trhitclus->GetAdcSum();
    event.sftlayerc.push_back(layer);
    event.sftclssize.push_back(clustersize);
    event.sftclssizelx.push_back(clusterlxsize);
    event.sftclssizelz.push_back(clusterlzsize);
    event.sftposlxc.push_back(lx);
    event.sftposgzc.push_back(gz);
    event.sftposgxc.push_back(globalx);
    event.sftposgyc.push_back(globaly);

    cluster.sftlayerc = layer;
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
      std::cout << __FILE__ << "  " << __LINE__ << " layer  " << layer << std::endl;
      std::cout << __FILE__ << "  " << __LINE__ << " lx  " << lx << std::endl;
      std::cout << __FILE__ << "  " << __LINE__ << " global z  " << gz << std::endl;
      std::cout << __FILE__ << "  " << __LINE__ << " global x  " << globalx << std::endl;
      std::cout << __FILE__ << "  " << __LINE__ << " global y  " << globaly << std::endl;
    }
  }


  if(is_doTracking){  
    TrAna->TrackSearchSFTT();//making index



    //This loop is for evaluating tracking performance while event matching between
    //primary info and reconstructed track info.
    int nt=TrAna->GetNtracksSFTT();
    event.ntr=nt;
    if(Verbosity>10){
      std::cout << __FILE__ << "  " << __LINE__ << " ntrack  " << nt << std::endl;
    }
    for( int it=0; it<nt; ++it ){
      TrLocalTrack *tp=TrAna->GetTrackSFTT(it);
      int nh=tp->GetNHit();
      double chisqr=tp->GetChiSquare();
      //double x0=tp->GetX0(), y0=tp->GetY0();
      double u0=tp->GetU0(), v0=tp->GetV0();//fitting results of u0 and v0 position
                                            //TODO : what is the definition of u0 and v0 ?
      double xtgt=tp->GetX( 0.0 ), ytgt=tp->GetY( 0.0 );
      double utgt=u0, vtgt=v0;
      
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
    TTree *ctree = dynamic_cast<TTree *>(gFile->Get("ctree"));
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
    TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
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
  int nbin = 4096;
  for(int i = 0; i < totalch; ++i) {
    adcHigh[i] = new TH1I(Form("ADC_HIGH_%d", i),
        Form("ADC high gain %d", i),
        nbin, 0, 4096);
    adcHigh[i]->SetXTitle("ch.");
    tdcLeading[i] = new TH1I(Form("TDC_LEADING_%d", i),
        Form("TDC leading %d", i),
        nbin, 0, 4096);
    tdcLeading[i]->SetXTitle("ch.");
    tdcTrailing[i] = new TH1I(Form("TDC_TRAILING_%d", i),
        Form("TDC trailing %d", i),
        nbin, 0, 4096);
    tdcTrailing[i]->SetXTitle("ch.");
  }

  for(int ich = 0; ich < totalch; ich++){
    pede_mean[ich] = pedestalmean_[ich];
    pede_sigma[ich] = pedestalsigma_[ich];
  }
  TH1F* hitprofile_ch[4][3];//layer , type(xuv)

  //cluster information
  TH1F* hitprofile[4][3];//layer , type(xuv)
  TH1I* hitmulti[4][3];//layer , type(xuv)

  const char XUV[3][2] = {"X","U","V"};
  for(int ilr =0 ; ilr<4; ilr++){
    for(int itype=0;itype<3;itype++){
      hitprofile[ilr][itype] = new TH1F(Form("Hitproflayer%s%d",XUV[itype],ilr),Form("Hit profile %s layer %d",XUV[itype],ilr), 40,-10,10); 
      hitprofile[ilr][itype]->SetXTitle("local pos. [mm]");
    }
  }
  
  for(int ilr =0 ; ilr<4; ilr++){
    for(int itype=0;itype<3;itype++){
      hitmulti[ilr][itype] = new TH1I(Form("HitMulti%s%d",XUV[itype],ilr),Form("Hit Multiplicity %s layer %d",XUV[itype],ilr), 20,0,20); 
      hitmulti[ilr][itype]->SetXTitle("hit multiplicity");
    }
  }



  std::cout << "Init "<< std::endl;
  return true;

}
