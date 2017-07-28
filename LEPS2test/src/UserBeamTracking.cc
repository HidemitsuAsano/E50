/*
  UserBeamTracking.cc

  2016/4 K.Shirotori 
*/

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
#include "PrimInfo.hh"
#include "SpecLib.hh"

#include <TStopwatch.h>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

//const int Verbosity = 100;
const int Verbosity = 0;

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventBeamTracking 
  : public VEvent
{

private:
  RawData *rawData;
  TrAnalyzer *TrAna;
  //  HodoAnalyzer *HodoAna;

public:
  EventBeamTracking();
  ~EventBeamTracking();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventBeamTracking::EventBeamTracking()
  : VEvent(),
    rawData(0),
    TrAna(new TrAnalyzer())
    //HodoAna(new HodoAnalyzer())
{
}

EventBeamTracking::~EventBeamTracking()
{
  //  if (HodoAna)   delete HodoAna;
  if (TrAna)   delete TrAna;
  if (rawData) delete rawData;
}

struct Event{
  unsigned int nevents;
  //Primary
  int    prinhits;
  double priposx, priposy, priposz;
  //pbeam : magnitude of beam mom. ubeam : beamMom.x()/beamMom.z() vbeam : beamMom.y()/beamMom.z() 
  double pbeam, ubeam, vbeam, abmom;
  double prim1, prip1;
  double pritheta1, priphi1, prithetacm1, priphicm1;
  double prim2, prip2;
  double pritheta2, priphi2, prithetacm2, priphicm2;
  double pxangle, pyangle;//true track angle

  //Generated Beam
  int    gbnhits;
  std::vector<int>    gbid, gbtype;
  std::vector<double> gbp, gbpx, gbpy, gbpz;
  std::vector<int>    gbpid;
  std::vector<double> gbvx, gbvy;

  //SFT (cluster based info.)
  int    sftnhits;
  std::vector<int> sftlayer;
  std::vector<double> sftposx, sftposy;//always 0 if realistic detector is choosed
  int    sftnclus;
  std::vector<int> sftlayerc;
  std::vector<int> sftclssize;
  std::vector<double> sftclssizelx, sftclssizelz;//
  std::vector<double> sftposlxc, sftposgzc;//
  std::vector<double> sftposgxc,sftposgyc;//global position calculated from localx and tiltangle
  std::vector<double> sftdl;//drift length ??? maybe, not necessary ?
  //T0
  int    t0nhits;
  std::vector<int>    t0layer, t0seg;
  std::vector<double> t0time, t0edep;
  std::vector<double> t0path, t0p;
  std::vector<double> t0posx, t0posy;
  std::vector<int>    t0pid;
  std::vector<double> t0beta;

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



bool EventBeamTracking::ProcessingBegin()
{
 return true;
}

bool EventBeamTracking::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->Decode(In) ) return false;
  //std::cout << "***" << std::endl;

  //**************************************************************************
  //******************RawData

  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  TTree *ctree = dynamic_cast<TTree *>(gFile->Get("ctree"));
  InitializeEvent();
  /*
  //PrimaryInfo
  {
    const PrimInfoContainer &cont=rawData->GetPrimHC();  
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      PrimInfo *hit=cont[i];
      event.priposx = hit->GetVertX();
      event.priposy = hit->GetVertY();
      event.priposz = hit->GetVertZ();
      
      event.pbeam = hit->GetBeamMom();
      event.ubeam = hit->GetBeamU();
      event.vbeam = hit->GetBeamV();
      event.pxangle= hit->GetBeamXangle();
      event.pyangle= hit->GetBeamYangle();
      event.abmom = hit->GetAnaBeamMom();
      
      event.prim1 = hit->GetMass1();
      event.prip1 = hit->GetMom1();
      event.pritheta1 = hit->GetTheta1();
      event.priphi1 = hit->GetPhi1();
      event.prithetacm1 = hit->GetThetaCM1();
      event.priphicm1 = hit->GetPhiCM1();
      
      event.prim2 = hit->GetMass2();
      event.prip2 = hit->GetMom2();
      event.pritheta2 = hit->GetTheta2();
      event.priphi2 = hit->GetPhi2();
      event.prithetacm2 = hit->GetThetaCM2();
      event.priphicm2 = hit->GetPhiCM2();
    }
  }*/ 

  //T0
  {
    const HodoRHitContainer &cont=rawData->GetT0RHC();
    int nh=cont.size();
    event.t0nhits = nh;
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.t0layer.push_back(layer);
      event.t0seg.push_back(seg);

      for( int j=0; j<nt; j++ ){
  	event.t0time.push_back(hit->GetTime(j));
  	event.t0edep.push_back(hit->GetEdep(j));
  	event.t0path.push_back(hit->GetPath(j));
  	event.t0p.push_back(hit->GetMom(j));
  	event.t0posx.push_back(hit->GetPosX(j));
  	event.t0posy.push_back(hit->GetPosY(j));
  	event.t0pid.push_back(hit->GetPid(j));
  	event.t0beta.push_back(hit->GetBeta(j));
      }
    }
  }

  //SFT
  {
    for( int layer=0; layer<NumOfLayersSFT; ++layer ){
      const SFTRawHitContainer &cont =rawData->GetSFTRawHitContainer(layer);
      int nh=cont.size();
      //std::cout << "N raw hit " << nh << std::endl;
      event.sftnhits = nh;
      for( int i=0; i<nh; ++i ){
  	SFTRawHit *hit=cont[i];
  	int nt = hit->GetSize();
	event.sftlayer.push_back(hit->LayerId());
  	for( int j=0; j<nt; j++ ) {
  	  event.sftlayer.push_back(hit->LayerId());
  	  event.sftposx.push_back(hit->GetPosX(j));//only stored if type 0 detector is selected
  	  event.sftposy.push_back(hit->GetPosY(j));//only stored if type 0 detector is selectedn
  	  event.sftdl.push_back(hit->GetDL(j));
  	}
      }
    }
  }


  //////////////////////////Tracking
  {  
    TrAna->DecodeSFTRawHits( rawData );
    TrAna->SortTrHits();//sort TrHits by segment ID
    //SFT
    //TODO check the TrHits are sorted or not
    TrAna->SFTClustering();//clustering, 
    //
    TrAna->TrackSearchSFTT();//making index
    
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
        ctree->Fill();
        
        if(Verbosity>3){
          std::cout << __FILE__ << "  " << __LINE__ << " event  " << event.nevents << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " cluster  " << iclus << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " layer  " << layer << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " lx  " << lx << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " global z  " << gz << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " globalx  " << globalx << std::endl;
          std::cout << __FILE__ << "  " << __LINE__ << " globaly  " << globaly << std::endl;
        }
      }
    }
    
    
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
      event.x0diff.push_back(xtgt - event.priposx);
      event.y0diff.push_back(ytgt - event.priposy);
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
          event.xanglediff.push_back(xangle-event.pxangle);
          event.yanglediff.push_back(xangle-event.pyangle);
        }
      }
    }
  }
  event.nevents++;
  tree->Fill();
  if(Verbosity>0){
    std::cout << __FILE__ << "  " << __LINE__ << " fill " << std::endl;
  }
  return true;
}

void EventBeamTracking::InitializeEvent( void )
{
  //PriInfo
  event.priposx = -9999.0;
  event.priposy = -9999.0;
  event.priposz = -9999.0;
  event.pbeam = -9999.0;
  event.ubeam = -9999.0;
  event.vbeam = -9999.0;
  event.abmom = -9999.0;
  event.prim1 = -9999.0;
  event.prip1 = -9999.0;
  event.pritheta1 = -9999.0;
  event.priphi1 = -9999.0;
  event.prithetacm1 = -9999.0;
  event.priphicm1 = -9999.0;
  event.prim2 = -9999.0;
  event.prip2 = -9999.0;
  event.pritheta2 = -9999.0;
  event.priphi2 = -9999.0;
  event.prithetacm2 = -9999.0;
  event.priphicm2 = -9999.0;
  event.pxangle = -9999.0;
  event.pyangle = -9999.0;

  //SFT
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
  event.sftdl.clear();
  event.sftxangle.clear();
  event.sftyangle.clear();

  //T0
  event.t0nhits = -1;
  event.t0layer.clear();
  event.t0seg.clear();
  event.t0time.clear();
  event.t0edep.clear();
  event.t0path.clear();
  event.t0p.clear();
  event.t0posx.clear();
  event.t0posy.clear();
  event.t0pid.clear();
  event.t0beta.clear();

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

bool EventBeamTracking::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventBeamTracking;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("ctree","cluster wise tree");
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
  
  tree->Branch("event", &event.nevents);
  //PriInfo  
  tree->Branch("priposx", &event.priposx);
  tree->Branch("priposy", &event.priposy);
  tree->Branch("priposz", &event.priposz);
  tree->Branch("pbeam", &event.pbeam);
  tree->Branch("ubeam", &event.ubeam);
  tree->Branch("vbeam", &event.vbeam);
  tree->Branch("abmom", &event.abmom);
  tree->Branch("prim1", &event.prim1);
  tree->Branch("prip1", &event.prip1);
  tree->Branch("pritheta1", &event.pritheta1);
  tree->Branch("priphi1", &event.priphi1);
  tree->Branch("prithetacm1", &event.prithetacm1);
  tree->Branch("priphicm1", &event.priphicm1);
  tree->Branch("prim2", &event.prim2);
  tree->Branch("prip2", &event.prip2);
  tree->Branch("pritheta2", &event.pritheta2);
  tree->Branch("priphi2", &event.priphi2);
  tree->Branch("prithetacm2", &event.prithetacm2);
  tree->Branch("priphicm2", &event.priphicm2);
  tree->Branch("pxangle",&event.pxangle);
  tree->Branch("pyangle",&event.pyangle);

  //SFT
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
  
  tree->Branch("sftdl", &event.sftdl);//drift length ???

  //T0
  tree->Branch("t0nhits", &event.t0nhits);
  tree->Branch("t0layer", &event.t0layer);
  tree->Branch("t0seg", &event.t0seg);
  tree->Branch("t0time", &event.t0time);
  tree->Branch("t0edep", &event.t0edep);
  tree->Branch("t0path", &event.t0path);
  tree->Branch("t0p", &event.t0p);
  tree->Branch("t0posx", &event.t0posx);
  tree->Branch("t0posy", &event.t0posy);
  tree->Branch("t0pid", &event.t0pid);
  tree->Branch("t0beta", &event.t0beta);
  
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

  return true;
}
