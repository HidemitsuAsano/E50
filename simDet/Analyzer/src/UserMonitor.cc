/*
  UserMonitor.cc

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
#include "PrimInfo.hh"
#include "HodoRawHit.hh"
#include "TrRawHit.hh"

 VEvent::VEvent()
 {
 }

 VEvent::~VEvent()
 {
 }

class EventMonitor 
  : public VEvent
{

private:
  RawData *rawData;

public:
  EventMonitor();
  ~EventMonitor();

  bool ProcessingBegin();
  bool ProcessingEnd();
  bool ProcessingNormal( std::ifstream & );
  bool InitializeHistograms(); 
  void InitializeEvent();
};

EventMonitor::EventMonitor()
  : VEvent(),
    rawData(0)
{
}

EventMonitor::~EventMonitor()
{
  if (rawData) delete rawData;
}

struct Event{
  //Primary
  int    prinhits;
  double priposx, priposy, priposz;
  double pbeam, ubeam, vbeam, abmom;
  double prim1, prip1;
  double pritheta1, priphi1, prithetacm1, priphicm1;
  double prim2, prip2;
  double pritheta2, priphi2, prithetacm2, priphicm2;

  //SFT
  int    sftnhits;
  std::vector<int>    sftlayer, sftseg;
  std::vector<double> sftposx, sftposy;
  std::vector<double> sftdl;

  //T0
  int    t0nhits;
  std::vector<int>    t0layer, t0seg;
  std::vector<double> t0time, t0edep;
  std::vector<double> t0path, t0p;
  std::vector<double> t0posx, t0posy;
  std::vector<int>    t0pid;
  std::vector<double> t0beta;

};
static Event event;

bool EventMonitor::ProcessingBegin()
{
 return true;
}

bool EventMonitor::ProcessingNormal( std::ifstream &In )
{
  const std::string funcname = "ProcessingNormal";

  rawData = new RawData;
  if( !rawData->DecodeRawHits(In) ) return false;
  //std::cout << "***" << std::endl;

  //**************************************************************************
  //******************RawData

  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

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
  }  

  //SFT
  {
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetSFTRawHitContainer(layer);
      int nh=cont.size();
      for( int i=0; i<nh; ++i ){
  	TrRawHit *hit=cont[i];
  	int nt = hit->GetSize();
  	event.sftnhits = nt;
  	for( int j=0; j<nt; j++ ) {
  	  event.sftlayer.push_back(hit->LayerId());
  	  event.sftseg.push_back(hit->WireId());
  	  event.sftposx.push_back(hit->GetPosX(j));
  	  event.sftposy.push_back(hit->GetPosY(j));
  	  event.sftdl.push_back(hit->GetDL(j));
  	}
      }
    }
  }

  //T0
  {
    const HodoRHitContainer &cont=rawData->GetT0RHC();
    int nh=cont.size();
    for( int i=0; i<nh; ++i ){
      HodoRawHit *hit=cont[i];
      int layer=hit->LayerId();
      int seg=hit->SegmentId();
      int nt = hit->GetSize();
      event.t0nhits = nt;
      for( int j=0; j<nt; j++ ){
  	event.t0layer.push_back(layer);
  	event.t0seg.push_back(seg);
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

  tree->Fill();

  return true;
}

void EventMonitor::InitializeEvent( void )
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

  //SFT
  event.sftnhits = -1;
  event.sftlayer.clear();
  event.sftseg.clear();
  event.sftposx.clear();
  event.sftposy.clear();
  event.sftdl.clear();

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
}


bool EventMonitor::ProcessingEnd()
{
  // gFile->Write();
  // gFile->Close();
  return true;
}

VEvent *ConfMan::EventAllocator()
{
  return new EventMonitor;
}

bool ConfMan:: InitializeHistograms()
{  
  HBTree("tree","tree of Spec");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));

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

  //SFT
  tree->Branch("sftnhits", &event.sftnhits);
  tree->Branch("sftlayer", &event.sftlayer);
  tree->Branch("sftseg", &event.sftseg);
  tree->Branch("sftposx", &event.sftposx);
  tree->Branch("sftposy", &event.sftposy);
  tree->Branch("sftdl", &event.sftdl);

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
  
  return true;
}
