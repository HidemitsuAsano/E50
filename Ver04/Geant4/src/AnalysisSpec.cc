/*
  AnalysisSpec.cc

  K.Shirotori
*/

#include "AnalysisSpec.hh"
#include "PrimaryInfo.hh"
#include "DataType.hh"

#include "TargetHit.hh"
#include "T0Hit.hh"
#include "BTrackerHit.hh"
#include "SITrackerHit.hh"
#include "AOTrackerHit.hh"
#include "IITrackerHit.hh"
#include "IOTrackerHit.hh"
#include "SOTrackerHit.hh"
#include "TofHit.hh"
#include "ITofHit.hh"
#include "PADHit.hh"
#include "RICHHit.hh"
#include "PID1Hit.hh"
#include "PID2Hit.hh"
#include "MuonFilterHit.hh"
#include "VDHit.hh"

#include "PidUtil.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include "Particle.hh"
#include "ParticleID.hh"
#include "RootHelper.hh"

#include <iomanip>
#include <vector>
#include <TTree.h>
#include "Random/Randomize.h"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

#define check 0

const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

AnalysisSpec::AnalysisSpec( const G4String & histname )
  : Analysis(histname)
{
  DefineHistograms();
}

AnalysisSpec::~AnalysisSpec()
{
}

void AnalysisSpec::BeginOfRun( const G4Run *aRun )
{
  G4cout << "AnalysisSpec::BeginOfRun" << G4endl;
  trigNum=0;
}

void AnalysisSpec::BeginOfEvent( const G4Event *anEvent )
{
  fTriggered = false;
}


void AnalysisSpec::EndOfRun( const G4Run *aRun )
{
  G4cout << "AnalysisSpec::EndOfRun::Total # of triggers ="
	 << std::setw(8) << trigNum << G4endl;
}

void AnalysisSpec::PrimaryGeneration( const PrimaryInfo *info )
{
  PrimaryInfo *pI = info->Copy();
  delete pInfo_; pInfo_=pI;
}

const G4double PionMass        = 0.1395701;
const G4double KaonMass        = 0.493677;
const G4double ProtonMass      = 0.93827200;
const G4double DMesonZMass     = 1.86480;
const G4double DMesonStarMMass = 2.010270;
const G4double LambdacPMass    = 2.286460;

#ifndef MaxHits 
#define MaxHits 30
#endif

struct Event{
  //Primary
  G4double priposx, priposy, priposz;
  G4double prim1, prim2;
  G4double prip1, prip2, prip1cm, prip2cm;
  G4double pritheta1, pritheta2, pritheta1cm, pritheta2cm;
  G4double priphi1, priphi2, priphi1cm, priphi2cm;
  G4double pricos1, pricos2, pricos1cm, pricos2cm;
  G4double priu1, priu2, priu1cm, priu2cm;
  G4double priv1, priv2, priv1cm, priv2cm;
  G4double pribeta1, pribeta2;
  G4double pbeam, ubeam, vbeam;
  G4double prip3, priu3, priv3;

  //Target
  G4int    tgtnhits;
  std::vector<G4int>    tgtpid, tgtid;
  std::vector<G4double> tgtedep, tgttime;
  std::vector<G4double> tgtp;
  std::vector<G4double> tgttheta, tgtphi;
  std::vector<G4double> tgtposx, tgtposy;
  std::vector<G4double> tgtpath, tgtbeta;
  std::vector<G4double> tgtmass;

  //T0
  G4int    t0nhits;
  std::vector<G4int>    t0layer, t0seg;
  std::vector<G4int>    t0pid, t0id;
  std::vector<G4double> t0edep, t0time;
  std::vector<G4double> t0p;
  std::vector<G4double> t0theta, t0phi;
  std::vector<G4double> t0posx, t0posy;
  std::vector<G4double> t0path, t0beta;
  std::vector<G4double> t0mass;

  //BTracker
  G4int    btrnhits;
  std::vector<G4int> btrlayer;
  std::vector<G4int>    btrpid, btrid;
  std::vector<G4double> btredep, btrtime;
  std::vector<G4double> btrp;
  std::vector<G4double> btrtheta, btrphi;
  std::vector<G4double> btrposx, btrposy;
  std::vector<G4double> btrpath, btrbeta;
  std::vector<G4double> btrmass;

  //SITracker
  G4int    sitrnhits;
  std::vector<G4int> sitrlayer;
  std::vector<G4int>    sitrpid, sitrid;
  std::vector<G4double> sitredep, sitrtime;
  std::vector<G4double> sitrp;
  std::vector<G4double> sitrtheta, sitrphi;
  std::vector<G4double> sitrposx, sitrposy;
  std::vector<G4double> sitrpath, sitrbeta;
  std::vector<G4double> sitrmass;

  //AOtracker
  G4int    aotrnhits;
  std::vector<G4int> aotrlayer;
  std::vector<G4int>    aotrpid, aotrid;
  std::vector<G4double> aotredep, aotrtime;
  std::vector<G4double> aotrp;
  std::vector<G4double> aotrtheta, aotrphi;
  std::vector<G4double> aotrposx, aotrposy;
  std::vector<G4double> aotrpath, aotrbeta;
  std::vector<G4double> aotrmass;

  //IITracker
  G4int    iitrnhits;
  std::vector<G4int> iitrlayer;
  std::vector<G4int>    iitrpid, iitrid;
  std::vector<G4double> iitredep, iitrtime;
  std::vector<G4double> iitrp;
  std::vector<G4double> iitrtheta, iitrphi;
  std::vector<G4double> iitrposx, iitrposy;
  std::vector<G4double> iitrpath, iitrbeta;
  std::vector<G4double> iitrmass;

  //IOTracker
  G4int    iotrnhits;
  std::vector<G4int> iotrlayer;
  std::vector<G4int>    iotrpid, iotrid;
  std::vector<G4double> iotredep, iotrtime;
  std::vector<G4double> iotrp;
  std::vector<G4double> iotrtheta, iotrphi;
  std::vector<G4double> iotrposx, iotrposy;
  std::vector<G4double> iotrpath, iotrbeta;
  std::vector<G4double> iotrmass;

  //SOTracker
  G4int    sotrnhits;
  std::vector<G4int> sotrlayer;
  std::vector<G4int>    sotrpid, sotrid;
  std::vector<G4double> sotredep, sotrtime;
  std::vector<G4double> sotrp;
  std::vector<G4double> sotrtheta, sotrphi;
  std::vector<G4double> sotrposx, sotrposy;
  std::vector<G4double> sotrpath, sotrbeta;
  std::vector<G4double> sotrmass;

  //Tof
  G4int    tofnhits;
  std::vector<G4int>    toflayer, tofseg;
  std::vector<G4int>    tofpid, tofid;
  std::vector<G4double> tofedep, toftime;
  std::vector<G4double> tofp;
  std::vector<G4double> toftheta, tofphi;
  std::vector<G4double> tofposx, tofposy;
  std::vector<G4double> tofpath, tofbeta;
  std::vector<G4double> tofmass;
  std::vector<G4double> tofm2_0, tofm2_1, tofm2_2;

  //ITof
  G4int    itofnhits;
  std::vector<G4int>    itoflayer, itofseg;
  std::vector<G4int>    itofpid, itofid;
  std::vector<G4double> itofedep, itoftime;
  std::vector<G4double> itofp;
  std::vector<G4double> itoftheta, itofphi;
  std::vector<G4double> itofposx, itofposy;
  std::vector<G4double> itofpath, itofbeta;
  std::vector<G4double> itofmass;
  std::vector<G4double> itofm2_0, itofm2_1, itofm2_2;

  //PAD
  G4int    padnhits;
  std::vector<G4int>    padlayer, padseg;
  std::vector<G4int>    padpid, padid;
  std::vector<G4double> padedep, padtime;
  std::vector<G4double> padp;
  std::vector<G4double> padtheta, padphi;
  std::vector<G4double> padposx, padposy;
  std::vector<G4double> padpath, padbeta;
  std::vector<G4double> padmass;
  std::vector<G4double> padm2_0, padm2_1, padm2_2;

  //RICH
  G4int    richnhits;
  std::vector<G4int>    richlayer, richseg;
  std::vector<G4int>    richpid, richid;
  std::vector<G4double> richedep, richtime;
  std::vector<G4double> richp;
  std::vector<G4double> richtheta, richphi;
  std::vector<G4double> richposx, richposy;
  std::vector<G4double> richpath, richbeta;
  std::vector<G4double> richmass;

  //PID1
  G4int    pid1nhits;
  std::vector<G4int>    pid1layer, pid1seg;
  std::vector<G4int>    pid1pid, pid1id;
  std::vector<G4double> pid1edep, pid1time;
  std::vector<G4double> pid1p;
  std::vector<G4double> pid1theta, pid1phi;
  std::vector<G4double> pid1posx, pid1posy;
  std::vector<G4double> pid1path, pid1beta;
  std::vector<G4double> pid1mass;

  //PID2
  G4int    pid2nhits;
  std::vector<G4int>    pid2layer, pid2seg;
  std::vector<G4int>    pid2pid, pid2id;
  std::vector<G4double> pid2edep, pid2time;
  std::vector<G4double> pid2p;
  std::vector<G4double> pid2theta, pid2phi;
  std::vector<G4double> pid2posx, pid2posy;
  std::vector<G4double> pid2path, pid2beta;
  std::vector<G4double> pid2mass;

  //MF
  G4int    mfnhits;
  std::vector<G4int>    mflayer, mfseg;
  std::vector<G4int>    mfpid, mfid;
  std::vector<G4double> mfedep, mftime;
  std::vector<G4double> mfp;
  std::vector<G4double> mftheta, mfphi;
  std::vector<G4double> mfposx, mfposy;
  std::vector<G4double> mfpath, mfbeta;
  std::vector<G4double> mfmass;

  //VD
  G4int    vdnhits;
  std::vector<G4int> vdlayer;
  std::vector<G4int>    vdpid, vdid;
  std::vector<G4double> vdedep, vdtime;
  std::vector<G4double> vdp;
  std::vector<G4double> vdtheta, vdphi;
  std::vector<G4double> vdposx, vdposy;
  std::vector<G4double> vdpath, vdbeta;
  std::vector<G4double> vdmass;

  //PID
  G4int ntKp,  ntKm;
  G4int ntPip, ntPim;
  G4int ntP,   ntPb;
  G4int ntMup, ntMum;
  G4int ntEp,  ntEm;

  std::vector<G4double> pkp,  ukp,  vkp,  thetakp,  phikp;  
  std::vector<G4double> pkm,  ukm,  vkm,  thetakm,  phikm;  
  std::vector<G4double> ppip, upip, vpip, thetapip, phipip;;
  std::vector<G4double> ppim, upim, vpim, thetapim, phipim;;
  std::vector<G4double> pp,   up,   vp,   thetap,   phip;   
  std::vector<G4double> ppb,  upb,  vpb,  thetapb,  phipb;  
  std::vector<G4double> pmup, umup, vmup, thetamup, phimup;
  std::vector<G4double> pmum, umum, vmum, thetamum, phimum;
  std::vector<G4double> pep,  uep,  vep,  thetaep,  phiep;  
  std::vector<G4double> pem,  uem,  vem,  thetaem,  phiem;

  //PID flag
  std::vector<G4int> ntdet;
  G4int ntrich,  ntpid1, ntpid2, nttof, ntitof, ntpad, ntmf;
  std::vector<G4int> fkp, fkm, fpip, fpim, fp, fpb, fmup, fmum, fep, fem;

  //Output
  G4double pkpx[MaxHits],  pkpy[MaxHits],  pkpz[MaxHits];
  G4double pkmx[MaxHits],  pkmy[MaxHits],  pkmz[MaxHits];
  G4double ppipx[MaxHits], ppipy[MaxHits], ppipz[MaxHits];
  G4double ppimx[MaxHits], ppimy[MaxHits], ppimz[MaxHits];
  G4double ppx[MaxHits],   ppy[MaxHits],   ppz[MaxHits];
  G4double ppbx[MaxHits],  ppby[MaxHits],  ppbz[MaxHits];
  G4double pmupx[MaxHits], pmupy[MaxHits], pmupz[MaxHits];
  G4double pmumx[MaxHits], pmumy[MaxHits], pmumz[MaxHits];
  G4double pepx[MaxHits],  pepy[MaxHits],  pepz[MaxHits];
  G4double pemx[MaxHits],  pemy[MaxHits],  pemz[MaxHits];

  G4int ntB;
  G4double pB[MaxHits],  uB[MaxHits],  vB[MaxHits];
  G4double pBx[MaxHits], pBy[MaxHits], pBz[MaxHits];

  G4int flkp[MaxHits], flkm[MaxHits];
  G4int flpip[MaxHits], flpim[MaxHits];
  G4int flp[MaxHits], flpb[MaxHits];
  G4int flmup[MaxHits], flmum[MaxHits];
  G4int flep[MaxHits], flem[MaxHits];
};
static Event event;

void AnalysisSpec::EndOfEvent( const G4Event *anEvent )
{
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
  InitializeEvent();

  // Primary Information
  std::vector <G4ThreeVector> BeamPCont, BeamXCont;
  G4ThreeVector primPos;
  if( fActive_ && pInfo_ ){
    G4double x=pInfo_->x/mm, y=pInfo_->y/mm, z=pInfo_->z/mm; 
    G4double m1=pInfo_->m1/GeV, m2=pInfo_->m2/GeV;
    G4double p1=pInfo_->p1/GeV, p2=pInfo_->p2/GeV;
    G4double p1cm=pInfo_->p1cm/GeV, p2cm=pInfo_->p2cm/GeV;
    G4double theta1=pInfo_->theta1/degree, theta2=pInfo_->theta2/degree;
    G4double theta1cm=pInfo_->theta1cm/degree, theta2cm=pInfo_->theta2cm/degree;
    G4double phi1=pInfo_->phi1/degree, phi2=pInfo_->phi2/degree;
    G4double phi1cm=pInfo_->phi1cm/degree, phi2cm=pInfo_->phi2cm/degree;
    G4double u1=pInfo_->u1, u2=pInfo_->u2;
    G4double u1cm=pInfo_->u1cm, u2cm=pInfo_->u2cm;
    G4double v1=pInfo_->v1, v2=pInfo_->v2;
    G4double v1cm=pInfo_->v1cm, v2cm=pInfo_->v2cm;
    G4double beta1=pInfo_->beta1, beta2=pInfo_->beta2;
    G4double pb=pInfo_->pbeam/GeV, ub=pInfo_->ubeam, vb=pInfo_->vbeam;
    G4double cos1=cos(pInfo_->theta1), cos2=cos(pInfo_->theta2);
    G4double cos1cm=cos(pInfo_->theta1cm), cos2cm=cos(pInfo_->theta2cm);
    G4double p3=pInfo_->p3/GeV, u3=pInfo_->u3, v3=pInfo_->v3;

    event.priposx = x, event.priposy = y, event.priposz = z; 
    event.prim1 = m1, event.prim2 = m2; 
    event.prip1 = p1, event.prip2 = p2; 
    event.prip1cm = p1cm, event.prip2cm = p2cm; 
    event.pritheta1 = theta1, event.pritheta2 = theta2; 
    event.pritheta1cm = theta1cm, event.pritheta2cm = theta2cm; 
    event.priphi1 = phi1, event.priphi2 = phi2; 
    event.priphi1cm = phi1cm, event.priphi2cm = phi2cm; 
    event.pricos1 = cos1, event.pricos2 = cos2; 
    event.pricos1cm = cos1cm, event.pricos2cm = cos2cm; 
    event.priu1 = u1, event.priu2 = u2; 
    event.priu1cm = u1cm, event.priu2cm = u2cm; 
    event.priv1 = v1, event.priv2 = v2; 
    event.priv1cm = v1cm, event.priv2cm = v2cm; 
    event.pribeta1 = beta1, event.pribeta2 = beta2; 
    event.pbeam = pb, event.ubeam = ub, event.vbeam = vb; 
    event.prip3 = p3, event.priu3 = u3, event.priv3 = v3; 

    G4double pt=pb/sqrt(1.+ub*ub+vb*vb);
    G4ThreeVector Pos( x, y, z );
    G4ThreeVector bMom( pt*ub, pt*vb, pt );
    primPos=Pos;
    BeamPCont.push_back(bMom); BeamXCont.push_back(primPos);
  }

#if 0
  PrintHitsInformation( anEvent, G4cout );
#endif

  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  G4int nhTgt=0;
  G4int nhT0=0, nhBTr=0, nhSITr=0, nhAOTr=0, nhIITr=0, nhIOTr=0, nhSOTr=0;
  G4int nhTof=0, nhITof=0, nhPAD=0, nhRICH=0, nhPID1=0, nhPID2=0, nhMF=0;
  G4int nhVD=0;
  TargetHitsCollection     *TgtHC=0;
  T0HitsCollection         *T0HC=0;
  BTrackerHitsCollection   *BTrHC=0;
  SITrackerHitsCollection  *SITrHC=0;
  AOTrackerHitsCollection  *AOTrHC=0;
  IITrackerHitsCollection  *IITrHC=0;
  IOTrackerHitsCollection  *IOTrHC=0;
  SOTrackerHitsCollection  *SOTrHC=0;
  TofHitsCollection        *TofHC=0;
  ITofHitsCollection       *ITofHC=0;
  PADHitsCollection        *PADHC=0;
  RICHHitsCollection       *RICHHC=0;
  PID1HitsCollection       *PID1HC=0;
  PID2HitsCollection       *PID2HC=0;
  MuonFilterHitsCollection *MFHC=0;
  VDHitsCollection         *VDHC=0;

  if(HCE){
    G4int colIdTgt=SDMan->GetCollectionID( "TargetCollection" );
    if(colIdTgt>=0){
      TgtHC=dynamic_cast<TargetHitsCollection *>( HCE->GetHC( colIdTgt ) );
      if(TgtHC) nhTgt=TgtHC->entries();
    }
    G4int colIdT0=SDMan->GetCollectionID( "T0Collection" );
    if(colIdT0>=0){
      T0HC=dynamic_cast<T0HitsCollection *>( HCE->GetHC( colIdT0 ) );
      if(T0HC) nhT0=T0HC->entries();
    }
    G4int colIdBTr=SDMan->GetCollectionID( "BTrackerCollection" );
    if(colIdBTr>=0){
      BTrHC=dynamic_cast<BTrackerHitsCollection *>( HCE->GetHC( colIdBTr ) );
      if(BTrHC) nhBTr=BTrHC->entries();
    }
    G4int colIdSITr=SDMan->GetCollectionID( "SITrackerCollection" );
    if(colIdSITr>=0){
      SITrHC=dynamic_cast<SITrackerHitsCollection *>( HCE->GetHC( colIdSITr ) );
      if(SITrHC) nhSITr=SITrHC->entries();
    }
    G4int colIdAOTr=SDMan->GetCollectionID( "AOTrackerCollection" );
    if(colIdAOTr>=0){
      AOTrHC=dynamic_cast<AOTrackerHitsCollection *>( HCE->GetHC( colIdAOTr ) );
      if(AOTrHC) nhAOTr=AOTrHC->entries();
    }
    G4int colIdIITr=SDMan->GetCollectionID( "IITrackerCollection" );
    if(colIdIITr>=0){
      IITrHC=dynamic_cast<IITrackerHitsCollection *>( HCE->GetHC( colIdIITr ) );
      if(IITrHC) nhIITr=IITrHC->entries();
    }
    G4int colIdIOTr=SDMan->GetCollectionID( "IOTrackerCollection" );
    if(colIdIOTr>=0){
      IOTrHC=dynamic_cast<IOTrackerHitsCollection *>( HCE->GetHC( colIdIOTr ) );
      if(IOTrHC) nhIOTr=IOTrHC->entries();
    }
    G4int colIdSOTr=SDMan->GetCollectionID( "SOTrackerCollection" );
    if(colIdSOTr>=0){
      SOTrHC=dynamic_cast<SOTrackerHitsCollection *>( HCE->GetHC( colIdSOTr ) );
      if(SOTrHC) nhSOTr=SOTrHC->entries();
    }
    G4int colIdTof=SDMan->GetCollectionID( "TofCollection" ); 
    if(colIdTof>=0){
      TofHC=dynamic_cast<TofHitsCollection *>( HCE->GetHC( colIdTof ) ); 
      if(TofHC) nhTof=TofHC->entries();
    }
    G4int colIdITof=SDMan->GetCollectionID( "ITofCollection" ); 
    if(colIdITof>=0){
      ITofHC=dynamic_cast<ITofHitsCollection *>( HCE->GetHC( colIdITof ) ); 
      if(ITofHC) nhITof=ITofHC->entries();
    }
    G4int colIdPAD=SDMan->GetCollectionID( "PADCollection" );
    if(colIdPAD>=0){
      PADHC=dynamic_cast<PADHitsCollection *>( HCE->GetHC( colIdPAD ) );
      if(PADHC) nhPAD=PADHC->entries();
    }
    G4int colIdRICH=SDMan->GetCollectionID( "RICHCollection" );
    if(colIdRICH>=0){
      RICHHC=dynamic_cast<RICHHitsCollection *>( HCE->GetHC( colIdRICH ) );
      if(RICHHC) nhRICH=RICHHC->entries();
    }
    G4int colIdPID1=SDMan->GetCollectionID( "PID1Collection" );
    if(colIdPID1>=0){
      PID1HC=dynamic_cast<PID1HitsCollection *>( HCE->GetHC( colIdPID1 ) );
      if(PID1HC) nhPID1=PID1HC->entries();
    }
    G4int colIdPID2=SDMan->GetCollectionID( "PID2Collection" );
    if(colIdPID2>=0){
      PID2HC=dynamic_cast<PID2HitsCollection *>( HCE->GetHC( colIdPID2 ) );
      if(PID2HC) nhPID2=PID2HC->entries();
    }
    G4int colIdMF=SDMan->GetCollectionID( "MuonFilterCollection" );
    if(colIdMF>=0){
      MFHC=dynamic_cast<MuonFilterHitsCollection *>( HCE->GetHC( colIdMF ) );
      if(MFHC) nhMF=MFHC->entries();
    }
    G4int colIdVD=SDMan->GetCollectionID( "VDCollection" );
    if(colIdVD>=0){
      VDHC=dynamic_cast<VDHitsCollection *>( HCE->GetHC( colIdVD ) );
      if(VDHC) nhVD=VDHC->entries();
    }
  }

  ConfMan *confMan = ConfMan::GetConfManager();
  G4bool counterHits = confMan->DetRootFill();

  if( counterHits ){
    ///////////////////////////////////////////////////////////////////
    ////////Counter Information
    ///////////////////////////

    //Target  
    event.tgtnhits = nhTgt;
    for( G4int i=0; i<nhTgt; i++ ){
      TargetHit *aHit=(*TgtHC)[i];

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.y(), mom.z());
      G4ThreeVector pos=aHit->GetPos();
      G4double xl=(-1.)*aHit->GetXLocal(), yl=(-1.)*aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //     G4cout<< i << " id=" << id
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;

      event.tgtpid.push_back(Pid);
      event.tgtid.push_back(id);
      event.tgtedep.push_back(edep/MeV);
      event.tgttime.push_back(time/ns);
      event.tgtp.push_back(p/GeV);
      event.tgttheta.push_back(theta);
      event.tgtphi.push_back(phi);
      event.tgtposx.push_back(xl/mm);
      event.tgtposy.push_back(yl/mm);
      event.tgtpath.push_back(path/mm);
      event.tgtbeta.push_back(beta);
      event.tgtmass.push_back(mass);
    }

    //T0
    event.t0nhits = nhT0;
    for( G4int i=0; i<nhT0; i++ ){
      T0Hit *aHit=(*T0HC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();

      //     G4cout<< i << " id=" << id << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;

      event.t0layer.push_back(layer);
      event.t0seg.push_back(seg);
      event.t0pid.push_back(Pid);
      event.t0id.push_back(id);
      event.t0edep.push_back(edep/MeV);
      event.t0time.push_back(time/ns);
      event.t0p.push_back(p/GeV);
      event.t0theta.push_back(theta);
      event.t0phi.push_back(phi);
      event.t0posx.push_back(xl/mm);
      event.t0posy.push_back(yl/mm);
      event.t0path.push_back(path/mm);
      event.t0beta.push_back(beta);
      event.t0mass.push_back(mass);
    }

    //BTracker
    event.btrnhits = nhBTr;
    for( G4int i=0; i<nhBTr; i++ ){
      BTrackerHit *aHit=(*BTrHC)[i];
      G4int layer=aHit->GetLayerID();
      
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //         G4cout<< i << " id=" << id << " Layer=" << layer
      //         	  << " Mass=" << mass 
      //     	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;
      
      event.btrlayer.push_back(layer);
      event.btrpid.push_back(Pid);
      event.btrid.push_back(id);
      event.btredep.push_back(edep/MeV);
      event.btrtime.push_back(time/ns);
      event.btrp.push_back(p/GeV);
      event.btrtheta.push_back(theta);
      event.btrphi.push_back(phi);
      event.btrposx.push_back(xl/mm);
      event.btrposy.push_back(yl/mm);
      event.btrpath.push_back(path/mm);
      event.btrbeta.push_back(beta);
      event.btrmass.push_back(mass);
    }
    
    //SITracker
    event.sitrnhits = nhSITr;
    for( G4int i=0; i<nhSITr; i++ ){
      SITrackerHit *aHit=(*SITrHC)[i];
      
      G4int layer=aHit->GetLayerID();
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //         G4cout<< i << " id=" << id << " Layer=" << layer
      //         	  << " Mass=" << mass 
      //     	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;
      
      event.sitrlayer.push_back(layer);
      event.sitrpid.push_back(Pid);
      event.sitrid.push_back(id);
      event.sitredep.push_back(edep/MeV);
      event.sitrtime.push_back(time/ns);
      event.sitrp.push_back(p/GeV);
      event.sitrtheta.push_back(theta);
      event.sitrphi.push_back(phi);
      event.sitrposx.push_back(xl/mm);
      event.sitrposy.push_back(yl/mm);
      event.sitrpath.push_back(path/mm);
      event.sitrbeta.push_back(beta);
      event.sitrmass.push_back(mass);
    }
    
    //AOTracker
    event.aotrnhits = nhAOTr;
    for( G4int i=0; i<nhAOTr; i++ ){
      AOTrackerHit *aHit=(*AOTrHC)[i];
      
      G4int layer=aHit->GetLayerID();
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //         G4cout<< i << " id=" << id << " Layer=" << layer
      //         	  << " Mass=" << mass 
      //     	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;
      
      event.aotrlayer.push_back(layer);
      event.aotrpid.push_back(Pid);
      event.aotrid.push_back(id);
      event.aotredep.push_back(edep/MeV);
      event.aotrtime.push_back(time/ns);
      event.aotrp.push_back(p/GeV);
      event.aotrtheta.push_back(theta);
      event.aotrphi.push_back(phi);
      event.aotrposx.push_back(xl/mm);
      event.aotrposy.push_back(yl/mm);
      event.aotrpath.push_back(path/mm);
      event.aotrbeta.push_back(beta);
      event.aotrmass.push_back(mass);
    }
    
    //IITracker 
    event.iitrnhits = nhIITr;
    for( G4int i=0; i<nhIITr; i++ ){
      IITrackerHit *aHit=(*IITrHC)[i];
      G4int layer=aHit->GetLayerID();
      
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      // G4cout<< i << " id=" << id << " Layer=" << layer
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;
      
      event.iitrlayer.push_back(layer);
      event.iitrpid.push_back(Pid);
      event.iitrid.push_back(id);
      event.iitredep.push_back(edep/MeV);
      event.iitrtime.push_back(time/ns);
      event.iitrp.push_back(p/GeV);
      event.iitrtheta.push_back(theta);
      event.iitrphi.push_back(phi);
      event.iitrposx.push_back(xl/mm);
      event.iitrposy.push_back(yl/mm);
      event.iitrpath.push_back(path/mm);
      event.iitrbeta.push_back(beta);
      event.iitrmass.push_back(mass);
    }

    //IOTracker 
    event.iotrnhits = nhIOTr;
    for( G4int i=0; i<nhIOTr; i++ ){
      IOTrackerHit *aHit=(*IOTrHC)[i];
      G4int layer=aHit->GetLayerID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      // G4cout<< i << " id=" << id << " Layer=" << layer
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;

      event.iotrlayer.push_back(layer);
      event.iotrpid.push_back(Pid);
      event.iotrid.push_back(id);
      event.iotredep.push_back(edep/MeV);
      event.iotrtime.push_back(time/ns);
      event.iotrp.push_back(p/GeV);
      event.iotrtheta.push_back(theta);
      event.iotrphi.push_back(phi);
      event.iotrposx.push_back(xl/mm);
      event.iotrposy.push_back(yl/mm);
      event.iotrpath.push_back(path/mm);
      event.iotrbeta.push_back(beta);
      event.iotrmass.push_back(mass);
    }
    
    //SOTracker
    event.sotrnhits = nhSOTr;
    for( G4int i=0; i<nhSOTr; i++ ){
      SOTrackerHit *aHit=(*SOTrHC)[i];
      G4int layer=aHit->GetLayerID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //         G4cout<< i << " id=" << id << " Layer=" << layer
      //         	  << " Mass=" << mass 
      //     	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;

      event.sotrlayer.push_back(layer);
      event.sotrpid.push_back(Pid);
      event.sotrid.push_back(id);
      event.sotredep.push_back(edep/MeV);
      event.sotrtime.push_back(time/ns);
      event.sotrp.push_back(p/GeV);
      event.sotrtheta.push_back(theta);
      event.sotrphi.push_back(phi);
      event.sotrposx.push_back(xl/mm);
      event.sotrposy.push_back(yl/mm);
      event.sotrpath.push_back(path/mm);
      event.sotrbeta.push_back(beta);
      event.sotrmass.push_back(mass);
    }
    
    //TOF  
    event.tofnhits = nhTof;
    for( G4int i=0; i<nhTof; i++ ){
      TofHit *aHit=(*TofHC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;
      
      event.toflayer.push_back(layer);
      event.tofseg.push_back(seg);
      event.tofpid.push_back(Pid);
      event.tofid.push_back(id);
      event.tofedep.push_back(edep/MeV);
      event.toftime.push_back(time/ns);
      event.tofp.push_back(p/GeV);
      event.toftheta.push_back(theta);
      event.tofphi.push_back(phi);
      event.tofposx.push_back(xl/mm);
      event.tofposy.push_back(yl/mm);
      event.tofpath.push_back(path/mm);
      event.tofbeta.push_back(beta);
      event.tofmass.push_back(mass);
      
      G4double C = 2.99792458E+8;
      G4double T0, V0, B0, m2_0;
      G4double T1, V1, B1, m2_1;
      G4double T2, V2, B2, m2_2;
      
      G4double P  = p*(1. + CLHEP::RandGauss::shoot( 0.0, 0.005 ) )/1000.;
      
      T0 = time + CLHEP::RandGauss::shoot( 0.0, 0.10 ); //nsec
      T1 = time + CLHEP::RandGauss::shoot( 0.0, 0.15 ); //nsec
      T2 = time + CLHEP::RandGauss::shoot( 0.0, 0.20 ); //nsec
      
      V0 = (path/1000.)/(T0*1.0E-9);
      V1 = (path/1000.)/(T1*1.0E-9);
      V2 = (path/1000.)/(T2*1.0E-9);
      
      B0 = V0/C;
      B1 = V1/C;
      B2 = V2/C;

      // std::cout << "B0= " << B0 << std::endl;
      // std::cout << "B1= " << B1 << std::endl;
      // std::cout << "B2= " << B2 << std::endl;

      m2_0 = (P/B0)*(P/B0)*(1.-B0*B0);
      m2_1 = (P/B1)*(P/B1)*(1.-B1*B1);
      m2_2 = (P/B2)*(P/B2)*(1.-B2*B2);
      
      event.tofm2_0.push_back(m2_0); 
      event.tofm2_1.push_back(m2_1); 
      event.tofm2_2.push_back(m2_2); 
    }
    
    //ITOF  
    event.itofnhits = nhITof;
    for( G4int i=0; i<nhITof; i++ ){
      ITofHit *aHit=(*ITofHC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();
      
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;
      
      event.itoflayer.push_back(layer);
      event.itofseg.push_back(seg);
      event.itofpid.push_back(Pid);
      event.itofid.push_back(id);
      event.itofedep.push_back(edep/MeV);
      event.itoftime.push_back(time/ns);
      event.itofp.push_back(p/GeV);
      event.itoftheta.push_back(theta);
      event.itofphi.push_back(phi);
      event.itofposx.push_back(xl/mm);
      event.itofposy.push_back(yl/mm);
      event.itofpath.push_back(path/mm);
      event.itofbeta.push_back(beta);
      event.itofmass.push_back(mass); 

      G4double C = 2.99792458E+8;
      G4double T0, V0, B0, m2_0;
      G4double T1, V1, B1, m2_1;
      G4double T2, V2, B2, m2_2;
      
      G4double P = p*(1. + CLHEP::RandGauss::shoot( 0.0, 0.01 ) )/1000.;
      
      T0 = time + CLHEP::RandGauss::shoot( 0.0, 0.07 ); //nsec
      T1 = time + CLHEP::RandGauss::shoot( 0.0, 0.10 ); //nsec
      T2 = time + CLHEP::RandGauss::shoot( 0.0, 0.15 ); //nsec
      
      V0 = (path/1000.)/(T0*1.0E-9);
      V1 = (path/1000.)/(T1*1.0E-9);
      V2 = (path/1000.)/(T2*1.0E-9);
      
      B0 = V0/C;
      B1 = V1/C;
      B2 = V2/C;
      
      // std::cout << "B0= " << B0 << std::endl;
      // std::cout << "B1= " << B1 << std::endl;
      // std::cout << "B2= " << B2 << std::endl;
      
      m2_0 = (P/B0)*(P/B0)*(1.-B0*B0);
      m2_1 = (P/B1)*(P/B1)*(1.-B1*B1);
      m2_2 = (P/B2)*(P/B2)*(1.-B2*B2);

      event.itofm2_0.push_back(m2_0); 
      event.itofm2_1.push_back(m2_1); 
      event.itofm2_2.push_back(m2_2); 
    }
    
    //PAD
    event.padnhits = nhPAD;
    for( G4int i=0; i<nhPAD; i++ ){
      PADHit *aHit=(*PADHC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();
      
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;

      event.padlayer.push_back(layer);
      event.padseg.push_back(seg);
      event.padpid.push_back(Pid);
      event.padid.push_back(id);
      event.padedep.push_back(edep/MeV);
      event.padtime.push_back(time/ns);
      event.padp.push_back(p/GeV);
      event.padtheta.push_back(theta);
      event.padphi.push_back(phi);
      event.padposx.push_back(xl/mm);
      event.padposy.push_back(yl/mm);
      event.padpath.push_back(path/mm);
      event.padbeta.push_back(beta);
      event.padmass.push_back(mass); 

      G4double C = 2.99792458E+8;
      G4double T0, V0, B0, m2_0;
      G4double T1, V1, B1, m2_1;
      G4double T2, V2, B2, m2_2;

      G4double P = p*(1. + CLHEP::RandGauss::shoot( 0.0, 0.01 ) )/1000.;
      
      T0 = time + CLHEP::RandGauss::shoot( 0.0, 0.07 ); //nsec
      T1 = time + CLHEP::RandGauss::shoot( 0.0, 0.10 ); //nsec
      T2 = time + CLHEP::RandGauss::shoot( 0.0, 0.15 ); //nsec
      
      V0 = (path/1000.)/(T0*1.0E-9);
      V1 = (path/1000.)/(T1*1.0E-9);
      V2 = (path/1000.)/(T2*1.0E-9);
      
      B0 = V0/C;
      B1 = V1/C;
      B2 = V2/C;

      // std::cout << "B0= " << B0 << std::endl;
      // std::cout << "B1= " << B1 << std::endl;
      // std::cout << "B2= " << B2 << std::endl;
      
      m2_0 = (P/B0)*(P/B0)*(1.-B0*B0);
      m2_1 = (P/B1)*(P/B1)*(1.-B1*B1);
      m2_2 = (P/B2)*(P/B2)*(1.-B2*B2);

      event.padm2_0.push_back(m2_0); 
      event.padm2_1.push_back(m2_1); 
      event.padm2_2.push_back(m2_2); 
    }

    //RICH
    event.richnhits = nhRICH;
    for( G4int i=0; i<nhRICH; i++ ){
      RICHHit *aHit=(*RICHHC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;
      
      event.richlayer.push_back(layer);
      event.richseg.push_back(seg);
      event.richpid.push_back(Pid);
      event.richid.push_back(id);
      event.richedep.push_back(edep/MeV);
      event.richtime.push_back(time/ns);
      event.richp.push_back(p/GeV);
      event.richtheta.push_back(theta);
      event.richphi.push_back(phi);
      event.richposx.push_back(xl/mm);
      event.richposy.push_back(yl/mm);
      event.richpath.push_back(path/mm);
      event.richbeta.push_back(beta);
      event.richmass.push_back(mass); 
    }
    
    //PID1
    event.pid1nhits = nhPID1;
    for( G4int i=0; i<nhPID1; i++ ){
      PID1Hit *aHit=(*PID1HC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;

      event.pid1layer.push_back(layer);
      event.pid1seg.push_back(seg);
      event.pid1pid.push_back(Pid);
      event.pid1id.push_back(id);
      event.pid1edep.push_back(edep/MeV);
      event.pid1time.push_back(time/ns);
      event.pid1p.push_back(p/GeV);
      event.pid1theta.push_back(theta);
      event.pid1phi.push_back(phi);
      event.pid1posx.push_back(xl/mm);
      event.pid1posy.push_back(yl/mm);
      event.pid1path.push_back(path/mm);
      event.pid1beta.push_back(beta);
      event.pid1mass.push_back(mass); 
    }

    //PID2
    event.pid2nhits = nhPID2;
    for( G4int i=0; i<nhPID2; i++ ){
      PID2Hit *aHit=(*PID2HC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;

      event.pid2layer.push_back(layer);
      event.pid2seg.push_back(seg);
      event.pid2pid.push_back(Pid);
      event.pid2id.push_back(id);
      event.pid2edep.push_back(edep/MeV);
      event.pid2time.push_back(time/ns);
      event.pid2p.push_back(p/GeV);
      event.pid2theta.push_back(theta);
      event.pid2phi.push_back(phi);
      event.pid2posx.push_back(xl/mm);
      event.pid2posy.push_back(yl/mm);
      event.pid2path.push_back(path/mm);
      event.pid2beta.push_back(beta);
      event.pid2mass.push_back(mass); 
    }
    
    //MF
    event.mfnhits = nhMF;
    for( G4int i=0; i<nhMF; i++ ){
      MuonFilterHit *aHit=(*MFHC)[i];
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();

      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
      
      //     G4cout<< i << " id=" << id 
      // 	  << " Layer=" << layer
      // 	  << " Seg=" << seg
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;
      
      event.mflayer.push_back(layer);
      event.mfseg.push_back(seg);
      event.mfpid.push_back(Pid);
      event.mfid.push_back(id);
      event.mfedep.push_back(edep/MeV);
      event.mftime.push_back(time/ns);
      event.mfp.push_back(p/GeV);
      event.mftheta.push_back(theta);
      event.mfphi.push_back(phi);
      event.mfposx.push_back(xl/mm);
      event.mfposy.push_back(yl/mm);
      event.mfpath.push_back(path/mm);
      event.mfbeta.push_back(beta);
      event.mfmass.push_back(mass); 
    }
    
    //VD  
    event.vdnhits = nhVD;
    for( G4int i=0; i<nhVD; i++ ){
      VDHit *aHit=(*VDHC)[i];
      G4int layer=aHit->GetLayerID();
      
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom(); 
      G4ThreeVector Mom(mom.x(), mom.z(), (-1.)*mom.y());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //     G4cout<< i << " id=" << id << " Layer=" << layer
      // 	  << " Mass=" << mass 
      // 	  << " PID=" << aHit->GetHitParticleName(0) << G4endl;
      
      event.vdlayer.push_back(layer);
      event.vdpid.push_back(Pid);
      event.vdid.push_back(id);
      event.vdedep.push_back(edep/MeV);
      event.vdtime.push_back(time/ns);
      event.vdp.push_back(p/GeV);
      event.vdtheta.push_back(theta);
      event.vdphi.push_back(phi);
      event.vdposx.push_back(xl/mm);
      event.vdposy.push_back(yl/mm);
      event.vdpath.push_back(path/mm);
      event.vdbeta.push_back(beta);
      event.vdmass.push_back(mass);
    }
  }

  G4bool trigCond1=false, trigCond2=false;
  //if( nhT0>=1 ){ trigCond1 = true;}
  //if( nhTof>=1 ){ trigCond1 = true;}
  //if( nhT0>=1 && nhTof>=1 && nhITof>=1 ){ trigCond1 = true; }
  //if( nhT0>=1 && nhTof>=2 && tofflag1 && nhITof>=1 ){ trigCond1 = true; }
  //if( nhT0>=1 && nhTof>=2 && tofflag1 && nhITof>=1 && stinflag1 && stinflag2 ){ trigCond1 = true; }
  
  trigCond1 = true;
  trigCond2 = true;

  if( DataFile_.is_open() && trigCond1 && trigCond2 )
    PrintHitsInformation( anEvent, DataFile_ );
  fTriggered = ( trigCond1 && trigCond2 );


  ///////////////////////////////////////////////////////////////////
  ////////Reaction analysis
  ///////////////////////////

  if( trigCond1 ) {
#if check
    std::cout<<"******************************"<<std::endl;
#endif

    //PID
    Particle *part = new Particle();

    G4int ntRich=0, ntPid1=0, ntPid2=0; 
    G4int ntTof=0, ntItof=0, ntPad=0, ntMF=0;

    G4int MinTrack    = 15;
    G4int MinTrackIn  = 10;
    G4int MinTrackOut =  7;
    G4int MinTrackPAD =  7;
    // G4int MinTrack    = 0;
    // G4int MinTrackIn  = 0;
    // G4int MinTrackOut = 0;
    // G4int MinTrackPAD = 0;

    for( G4int i=0; i<nhTgt; i++ ){
      TargetHit *aHittgt=(*TgtHC)[i];
      G4int id=aHittgt->GetTrackNo();
      G4int Pid=PIDFunc(aHittgt->GetHitParticleName(0));
      G4ThreeVector mom=aHittgt->GetMom();
      G4double p=mom.mag();    

      //TOF hit analysis
      for( G4int j=0; j<nhTof; j++ ){
	TofHit *aHit=(*TofHC)[j];
	G4int trackNo=aHit->GetTrackNo();
	
	if( id==trackNo ){
	  G4bool flag_track =false;
	  G4int trackhitsIn=0;
	  G4int trackhitsOut=0;
	  G4bool flag_rich =false;
	  G4bool flag_pid1 =false;

	  //Track hit In: SFT + IT1
	  for( G4int k=0; k<nhSITr; k++ ){
	    SITrackerHit *aHit1=(*SITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  for( G4int k=0; k<nhIITr; k++ ){
	    IITrackerHit *aHit1=(*IITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  //Track hit Out: ST1 + ST2
	  for( G4int k=0; k<nhSOTr; k++ ){
	    SOTrackerHit *aHit1=(*SOTrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsOut++; 
	  }
	  if( (trackhitsIn>=MinTrackIn && trackhitsOut>=MinTrackOut) 
	      || (trackhitsIn+trackhitsOut>=MinTrack) ) flag_track=true;
	  
	  //RICH ID
	  for( G4int k=0; k<nhRICH; k++ ){
	    RICHHit *aHit1=(*RICHHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_rich =true;
	  }
	  //PID1 ID
	  for( G4int k=0; k<nhPID1; k++ ){
	    PID1Hit *aHit1=(*PID1HC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid1 =true;
	  }

	  if( flag_track && flag_rich ){
#if check
	    std::cout<<"******1: RICH"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(1);
	    event.ntdet.push_back(1);
	    ntRich++;
	  }

	  if( flag_track && flag_pid1 && !flag_rich ){
#if check
	    std::cout<<"******2: PID1"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(2);
	    event.ntdet.push_back(2);
	    ntPid1++;
	  }

	  if( flag_track && !flag_pid1 && !flag_rich ){
#if check
	    std::cout<<"******3: TOF"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(3);
	    event.ntdet.push_back(3);
	    ntTof++;
	  }
	}
      }

      //ITOF hit analysis
      for( G4int j=0; j<nhITof; j++ ){
	ITofHit *aHit=(*ITofHC)[j];
	G4int trackNo=aHit->GetTrackNo();
	G4int layer=aHit->GetLayerID();
	
	if( id==trackNo ){
	  G4bool flag_track =false;
	  G4bool flag_track2 =false;
	  G4bool flag_aft =false;
	  G4int trackhitsIn=0;
	  G4int trackhitsOut=0;
	  G4bool flag_tof =false;
	  G4bool flag_pid1 =false;
	  G4bool flag_pid2 =false;

	  //Track hit In: SFT + IT1
	  for( G4int k=0; k<nhSITr; k++ ){
	    SITrackerHit *aHit1=(*SITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  for( G4int k=0; k<nhIITr; k++ ){
	    IITrackerHit *aHit1=(*IITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  //Track hit Out: ST1&2 or IT2 or AFT
	  for( G4int k=0; k<nhSOTr; k++ ){
	    SOTrackerHit *aHit1=(*SOTrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsOut++; 
	  }
	  for( G4int k=0; k<nhIOTr; k++ ){
	    IOTrackerHit *aHit1=(*IOTrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsOut++; 
	  }
	  for( G4int k=0; k<nhAOTr; k++ ){
	    AOTrackerHit *aHit1=(*AOTrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_aft=true; 
	  }
	  if( (trackhitsIn>=MinTrackIn && trackhitsOut>=MinTrackOut) 
	      || (trackhitsIn+trackhitsOut>=MinTrack)
	      || (trackhitsIn>=MinTrackIn && flag_aft) ) flag_track=true;
	  if( trackhitsIn>=MinTrackIn ) flag_track2=true;
	  
	  //TOF ID
	  for( G4int k=0; k<nhTof; k++ ){
	    TofHit *aHit1=(*TofHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_tof =true;
	  }
	  //PID1 ID
	  for( G4int k=0; k<nhPID1; k++ ){
	    PID1Hit *aHit1=(*PID1HC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid1 =true;
	  }
	  //PID2 ID
	  for( G4int k=0; k<nhPID2; k++ ){
	    PID2Hit *aHit1=(*PID2HC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid2 =true;
	  }

	  if( flag_track && !flag_tof && flag_pid1 ){
#if check
	    std::cout<<"******2: PID1"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(2);
	    event.ntdet.push_back(2);
	    ntPid1++;
	  }


	  if( flag_track && !flag_tof && !flag_pid1 && flag_pid2 ){
#if check
	    std::cout<<"******4: PID2"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(4);
	    event.ntdet.push_back(4);
	    ntPid2++;
	  }

	  if( (flag_track || flag_track2 ) && !flag_tof && !flag_pid1 && !flag_pid2 ){
#if check
	    std::cout<<"******5: ITOF"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     << " TrackOut= " << trackhitsOut
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(5);
	    event.ntdet.push_back(5);
	    ntItof++;
	  }
	}
      }

      //PAD hit analysis
      for( G4int j=0; j<nhPAD; j++ ){
	PADHit *aHit=(*PADHC)[j];
	G4int trackNo=aHit->GetTrackNo();
	
	if( id==trackNo ){
	  G4bool flag_track =false;
	  G4int trackhitsIn=0;
	  G4bool flag_tof =false;
	  G4bool flag_itof =false;
	  G4bool flag_rich =false;
	  G4bool flag_pid1 =false;
	  G4bool flag_pid2 =false;

	  //Track hit In: SFT + IT1
	  for( G4int k=0; k<nhSITr; k++ ){
	    SITrackerHit *aHit1=(*SITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  for( G4int k=0; k<nhIITr; k++ ){
	    IITrackerHit *aHit1=(*IITrHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) trackhitsIn++; 
	  }
	  if( trackhitsIn>=MinTrackPAD ) flag_track=true;

	  //TOF ID
	  for( G4int k=0; k<nhTof; k++ ){
	    TofHit *aHit1=(*TofHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_tof =true;
	  }
	  //ITOF ID
	  for( G4int k=0; k<nhITof; k++ ){
	    ITofHit *aHit1=(*ITofHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_itof =true;
	  }
	  //RICH ID
	  for( G4int k=0; k<nhRICH; k++ ){
	    RICHHit *aHit1=(*RICHHC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_rich =true;
	  }
	  //PID1 ID
	  for( G4int k=0; k<nhPID1; k++ ){
	    PID1Hit *aHit1=(*PID1HC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid1 =true;
	  }
	  //PID2 ID
	  for( G4int k=0; k<nhPID2; k++ ){
	    PID2Hit *aHit1=(*PID2HC)[k];
	    G4int trackNo1=aHit1->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid2 =true;
	  }
	  
	  if( flag_track && !flag_tof && !flag_itof && !flag_rich && !flag_pid1 && !flag_pid2 ){
#if check
	    std::cout<<"******6: PAD"<<std::endl;
	    std::cout<< "id= " << trackNo 
	    	     << " PID= " << Pid
	    	     << " Mom= " << p 
	    	     << " TrackIn= " << trackhitsIn
	    	     <<std::endl;
#endif
	    part->SetId(id);
	    part->SetPid(Pid);
	    part->SetMomX((-1.)*mom.y());
	    part->SetMomY(mom.z());
	    part->SetMomZ((-1.)*mom.x());
	    part->SetDetFlag(6);
	    event.ntdet.push_back(6);
	    ntPad++;
	  }
	}
      }

      //MuonFilter ID
      for( G4int j=0; j<nhMF; j++ ){
	MuonFilterHit *aHit=(*MFHC)[j];
	G4int trackNo=aHit->GetTrackNo(); 
	G4int layer=aHit->GetLayerID();
	G4int seg=aHit->GetSegmentID();
	G4double time=aHit->GetTime();

	if( id==trackNo && layer==1 ){
	  G4bool flag_rich=false;
	  G4bool flag_pid1=false;
	  G4bool flag_tof=false;
	  G4bool flag_itof=false;

	  for( G4int k=0; k<nhRICH; k++ ){
	    RICHHit *aHit=(*RICHHC)[k];
	    G4int trackNo1=aHit->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_rich=true;
	  }
	  for( G4int k=0; k<nhPID1; k++ ){
	    PID1Hit *aHit=(*PID1HC)[k];
	    G4int trackNo1=aHit->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_pid1=true;
	  }
	  for( G4int k=0; k<nhTof; k++ ){
	    TofHit *aHit=(*TofHC)[k];
	    G4int trackNo1=aHit->GetTrackNo();
	    if( trackNo==trackNo1 ) flag_tof=true;
	  }
	  if( (flag_rich || flag_pid1) && flag_tof ){ 
	    for( G4int l=0; l<nhMF; l++ ){
	      MuonFilterHit *aHit=(*MFHC)[l];
	      G4int trackNo2=aHit->GetTrackNo(); 
	      G4int layer2=aHit->GetLayerID();
	      G4int seg2=aHit->GetSegmentID();
	      G4double time2=aHit->GetTime();

	      if( trackNo==trackNo2 && layer2==2 ){
#if check		
		// std::cout<<"******7"<<std::endl;
		// std::cout<< "id= " << trackNo 
		// 	     << " PID= " << Pid
		// 	     << " Mom= " << p 
		// 	     << " Time= " << time 
		// 	     << " Seg= " << seg 
		// 	     << " Layer= " << layer 
		// 	     <<std::endl;
#endif
		part->SetId(id);
		part->SetPid(Pid);
		part->SetMomX((-1.)*mom.y());
		part->SetMomY(mom.z());
		part->SetMomZ((-1.)*mom.x());
		part->SetDetFlag(7);
		event.ntdet.push_back(7);
		ntMF++;
	      }
	    }
	  }
	}
      }
      /////////////////////
    } //PID end
  
    event.ntrich=ntRich;
    event.ntpid1=ntPid1;
    event.ntpid2=ntPid2;
    event.nttof =ntTof;
    event.ntitof=ntItof;
    event.ntpad =ntPad;
    event.ntmf  =ntMF;

    std::vector <G4ThreeVector> KaonPCont, KaonPXCont; 
    std::vector <G4ThreeVector> KaonMCont, KaonMXCont; 
    std::vector <G4ThreeVector> PionPCont, PionPXCont;
    std::vector <G4ThreeVector> PionMCont, PionMXCont;
    std::vector <G4ThreeVector> ProtonCont, ProtonXCont;
    std::vector <G4ThreeVector> ProtonbCont, ProtonbXCont;
    std::vector <G4ThreeVector> MuonPCont, MuonPXCont; 
    std::vector <G4ThreeVector> MuonMCont, MuonMXCont; 
    std::vector <G4ThreeVector> PositronCont, PositronXCont; 
    std::vector <G4ThreeVector> ElectronCont, ElectronXCont; 

    G4int ikp=0, ikm=0;
    G4int ipip=0, ipim=0;
    G4int ip=0, ipb=0;
    G4int iep=0, iem=0;
    G4int imup=0, imum=0;

    G4int ntPart = part->GetSize();

    for( G4int i=0; i<ntPart; i++ ){
      G4int pid=part->GetPid(i);
      G4double px=part->GetMomX(i), py=part->GetMomY(i), pz=part->GetMomZ(i);
      G4ThreeVector mom(px, py, pz);
      G4double p=mom.mag(), u=mom.x()/p, v=mom.y()/p;
      G4double theta = mom.theta()*Rad2Deg, phi=mom.phi()*Rad2Deg;
      G4int fldet = part->GetDetFlag(i);
	  
      if( pid==4 ){//K+
    	event.pkp.push_back(p/GeV);
    	event.ukp.push_back(u);
    	event.vkp.push_back(v);
    	event.thetakp.push_back(theta);	    
    	event.phikp.push_back(phi);
    	event.pkpx[ikp]=mom.x()/GeV;
    	event.pkpy[ikp]=mom.y()/GeV;
    	event.pkpz[ikp]=mom.z()/GeV;

	event.fkp.push_back(fldet);
	event.flkp[ikp]=fldet;
	
    	KaonPCont.push_back(mom);
	KaonPXCont.push_back(primPos);
    	ikp++;
      }
      if( pid==3 ){//K-
    	event.pkm.push_back(p/GeV);
    	event.ukm.push_back(u);
    	event.vkm.push_back(v);
    	event.thetakm.push_back(theta);	    
    	event.phikm.push_back(phi);
    	event.pkmx[ikm]=mom.x()/GeV;
    	event.pkmy[ikm]=mom.y()/GeV;
    	event.pkmz[ikm]=mom.z()/GeV;

	event.fkm.push_back(fldet);
	event.flkm[ikm]=fldet;
	      
    	KaonMCont.push_back(mom);
	KaonMXCont.push_back(primPos);
    	ikm++;
      }
      if( pid==2 ){//Pi+
    	event.ppip.push_back(p/GeV);
    	event.upip.push_back(u);
    	event.vpip.push_back(v);
    	event.thetapip.push_back(theta);	    
    	event.phipip.push_back(phi);
    	event.ppipx[ipip]=mom.x()/GeV;
    	event.ppipy[ipip]=mom.y()/GeV;
    	event.ppipz[ipip]=mom.z()/GeV;

	event.fpip.push_back(fldet);
	event.flpip[ipip]=fldet;
	
    	PionPCont.push_back(mom);
	PionPXCont.push_back(primPos);
    	ipip++;
      }
      if( pid==1 ){//Pi-
    	event.ppim.push_back(p/GeV);
    	event.upim.push_back(u);
    	event.vpim.push_back(v);
    	event.thetapim.push_back(theta);	    
    	event.phipim.push_back(phi);
    	event.ppimx[ipim]=mom.x()/GeV;
    	event.ppimy[ipim]=mom.y()/GeV;
    	event.ppimz[ipim]=mom.z()/GeV;

	event.fpim.push_back(fldet);
	event.flpim[ipim]=fldet;
	
    	PionMCont.push_back(mom);
	PionMXCont.push_back(primPos);
    	ipim++;
      }
      if( pid==10 ){//Proton
    	event.pp.push_back(p/GeV);
    	event.up.push_back(u);
    	event.vp.push_back(v);
    	event.thetap.push_back(theta);	    
    	event.phip.push_back(phi);
    	event.ppx[ip]=mom.x()/GeV;
    	event.ppy[ip]=mom.y()/GeV;
    	event.ppz[ip]=mom.z()/GeV;

	event.fp.push_back(fldet);
	event.flp[ip]=fldet;

    	ProtonCont.push_back(mom);
	ProtonXCont.push_back(primPos);
    	ip++;
      }
      if( pid==14 ){//Anti-Proton
    	event.ppb.push_back(p/GeV);
    	event.upb.push_back(u);
    	event.vpb.push_back(v);
    	event.thetapb.push_back(theta);	    
    	event.phipb.push_back(phi);
    	event.ppbx[ipb]=mom.x()/GeV;
    	event.ppby[ipb]=mom.y()/GeV;
    	event.ppbz[ipb]=mom.z()/GeV;

	event.fpb.push_back(fldet);
	event.flpb[ipb]=fldet;
	
    	ProtonbCont.push_back(mom);
	ProtonbXCont.push_back(primPos);
    	ipb++;
      }
      if( pid==6 ){ //Mu+
    	event.pmup.push_back(p/GeV);
    	event.umup.push_back(u);
    	event.vmup.push_back(v);
    	event.thetamup.push_back(theta);	    
    	event.phimup.push_back(phi);
    	event.pmupx[imup]=mom.x()/GeV;
    	event.pmupy[imup]=mom.y()/GeV;
    	event.pmupz[imup]=mom.z()/GeV;
	
	event.fmup.push_back(fldet);
	event.flmup[imup]=fldet;

    	MuonPCont.push_back(mom);
	MuonPXCont.push_back(primPos);
    	imup++;
      }
      if( pid==5 ){ //Mu-
    	event.pmum.push_back(p/GeV);
    	event.umum.push_back(u);
    	event.vmum.push_back(v);
    	event.thetamum.push_back(theta);	    
    	event.phimum.push_back(phi);
    	event.pmumx[imum]=mom.x()/GeV;
    	event.pmumy[imum]=mom.y()/GeV;
    	event.pmumz[imum]=mom.z()/GeV;

	event.fmum.push_back(fldet);
	event.flmum[imum]=fldet;
	
    	MuonMCont.push_back(mom);
	MuonMXCont.push_back(primPos);
    	imum++;
      }
      if( pid==8 ){ //e+
    	event.pep.push_back(p/GeV);
    	event.uep.push_back(u);
    	event.vep.push_back(v);
    	event.thetaep.push_back(theta);	    
    	event.phiep.push_back(phi);
    	event.pepx[iep]=mom.x()/GeV;
    	event.pepy[iep]=mom.y()/GeV;
    	event.pepz[iep]=mom.z()/GeV;

	event.fep.push_back(fldet);
	event.flep[iep]=fldet;
	
    	PositronCont.push_back(mom);
	PositronXCont.push_back(primPos);
    	iep++;
      }
      if( pid==7 ){ //e-
    	event.pem.push_back(p/GeV);
    	event.uem.push_back(u);
    	event.vem.push_back(v);
    	event.thetaem.push_back(theta);	    
    	event.phiem.push_back(phi);
    	event.pemx[iem]=mom.x()/GeV;
    	event.pemy[iem]=mom.y()/GeV;
    	event.pemz[iem]=mom.z()/GeV;

	event.fem.push_back(fldet);
	event.flem[iem]=fldet;
	
    	ElectronCont.push_back(mom);
	ElectronXCont.push_back(primPos);
    	iem++;
      }
    }

    //Number of Particles
    G4int nKp=KaonPCont.size();
    G4int nKm=KaonMCont.size();
    G4int nPip=PionPCont.size();  
    G4int nPim=PionMCont.size();
    G4int nP=ProtonCont.size();
    G4int nPb=ProtonbCont.size();
    G4int nMup=MuonPCont.size();  
    G4int nMum=MuonMCont.size();
    G4int nEp=PositronCont.size();  
    G4int nEm=ElectronCont.size();

    event.ntKp = nKp; 
    event.ntKm = nKm; 
    event.ntPip = nPip; 
    event.ntPim = nPim; 
    event.ntP = nP; 
    event.ntPb = nPb; 
    event.ntMup = nMup; 
    event.ntMum = nMum; 
    event.ntEp = nEp; 
    event.ntEm = nEm; 

    //Beam
    G4int ibe=0;
    G4int nB=BeamPCont.size();
    event.ntB = nB;
    for( G4int ib=0; ib<nB; ++ib ){
      G4ThreeVector pb=BeamPCont[ib];
      G4double ub=pb.x()/pb.z(), vb=pb.y()/pb.z();
      event.pB[ibe]=pb.mag();
      event.uB[ibe]=ub;
      event.vB[ibe]=vb;
      event.pBx[ibe]=pb.x();
      event.pBy[ibe]=pb.y();
      event.pBz[ibe]=pb.z();
      ibe++;
    }

    tree->Fill();  

    delete part;

    // //For RICH study
    // //G4cout<< "**************" << G4endl;
    // for( G4int i=0; i<nhVD; i++ ){
    //   VDHit *aHit=(*VDHC)[i]; 
    //   G4int layer=aHit->GetLayerID();
    //   G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
    //   G4ThreeVector mom=aHit->GetMom(); 
    //   G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();

    //   if( layer==2 ){
    // 	if( Pid==4 ){
    // 	  G4cout<< 321 <<G4endl;
    // 	  G4cout<< (-1.)*mom.y()/GeV << " " << mom.z()/GeV << " " << (-1.)*mom.x()/GeV <<G4endl;
    // 	  G4cout<< xl/mm << " " << yl/mm <<G4endl;
    // 	}
    // 	if( Pid==3 ){
    // 	  G4cout<< -321 <<G4endl;
    // 	  G4cout<< (-1.)*mom.y()/GeV << " " << mom.z()/GeV << " " << (-1.)*mom.x()/GeV <<G4endl;
    // 	  G4cout<< xl/mm << " " << yl/mm <<G4endl;
    // 	}
    // 	if( Pid==2 ){
    // 	  G4cout<< 211 <<G4endl;
    // 	  G4cout<< (-1.)*mom.y()/GeV << " " << mom.z()/GeV << " " << (-1.)*mom.x()/GeV <<G4endl;
    // 	  G4cout<< xl/mm << " " << yl/mm <<G4endl;
    // 	}
    // 	//if( Pid==1 && MaxMomPis<mom.mag()/GeV  ){
    // 	if( Pid==1 ){
    // 	  G4cout<< -211 <<G4endl;
    // 	  G4cout<< (-1.)*mom.y()/GeV << " " << mom.z()/GeV << " " << (-1.)*mom.x()/GeV <<G4endl;
    // 	  G4cout<< xl/mm << " " << yl/mm <<G4endl;
    // 	}
    // 	if( Pid==10 ){
    // 	  G4cout<< 2212 <<G4endl;
    // 	  G4cout<< (-1.)*mom.y()/GeV << " " << mom.z()/GeV << " " << (-1.)*mom.x()/GeV <<G4endl;
    // 	  G4cout<< xl/mm << " " << yl/mm <<G4endl;
    // 	}
    //   }
    // }
    // G4cout<< -1 << G4endl;
  }

  //////////////////////End
}

void AnalysisSpec::PrintHitsInformation( const G4Event *anEvent, std::ostream &ost ) const
{
  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  G4int nhTgt=0;
  G4int nhT0=0, nhBTr=0, nhSITr=0, nhAOTr=0, nhIITr=0, nhIOTr=0, nhSOTr=0;
  G4int nhTof=0, nhITof=0, nhPAD=0, nhRICH=0, nhPID1=0, nhPID2=0, nhMF=0;
  G4int nhVD=0;
  TargetHitsCollection     *TgtHC=0;
  T0HitsCollection         *T0HC=0;
  BTrackerHitsCollection   *BTrHC=0;
  SITrackerHitsCollection  *SITrHC=0;
  AOTrackerHitsCollection  *AOTrHC=0;
  IITrackerHitsCollection  *IITrHC=0;
  IOTrackerHitsCollection  *IOTrHC=0;
  SOTrackerHitsCollection  *SOTrHC=0;
  TofHitsCollection        *TofHC=0;
  ITofHitsCollection       *ITofHC=0;
  PADHitsCollection        *PADHC=0;
  RICHHitsCollection       *RICHHC=0;
  PID1HitsCollection       *PID1HC=0;
  PID2HitsCollection       *PID2HC=0;
  MuonFilterHitsCollection *MFHC=0;
  VDHitsCollection         *VDHC=0;

  if(HCE){
    G4int colIdTgt=SDMan->GetCollectionID( "TargetCollection" );
    if(colIdTgt>=0){
      TgtHC=dynamic_cast<TargetHitsCollection *>( HCE->GetHC( colIdTgt ) );
      if(TgtHC) nhTgt=TgtHC->entries();
    }
    G4int colIdT0=SDMan->GetCollectionID( "T0Collection" );
    if(colIdT0>=0){
      T0HC=dynamic_cast<T0HitsCollection *>( HCE->GetHC( colIdT0 ) );
      if(T0HC) nhT0=T0HC->entries();
    }
    G4int colIdBTr=SDMan->GetCollectionID( "BTrackerCollection" );
    if(colIdBTr>=0){
      BTrHC=dynamic_cast<BTrackerHitsCollection *>( HCE->GetHC( colIdBTr ) );
      if(BTrHC) nhBTr=BTrHC->entries();
    }
    G4int colIdSITr=SDMan->GetCollectionID( "SITrackerCollection" );
    if(colIdSITr>=0){
      SITrHC=dynamic_cast<SITrackerHitsCollection *>( HCE->GetHC( colIdSITr ) );
      if(SITrHC) nhSITr=SITrHC->entries();
    }
    G4int colIdAOTr=SDMan->GetCollectionID( "AOTrackerCollection" );
    if(colIdAOTr>=0){
      AOTrHC=dynamic_cast<AOTrackerHitsCollection *>( HCE->GetHC( colIdAOTr ) );
      if(AOTrHC) nhAOTr=AOTrHC->entries();
    }
    G4int colIdIITr=SDMan->GetCollectionID( "IITrackerCollection" );
    if(colIdIITr>=0){
      IITrHC=dynamic_cast<IITrackerHitsCollection *>( HCE->GetHC( colIdIITr ) );
      if(IITrHC) nhIITr=IITrHC->entries();
    }
    G4int colIdIOTr=SDMan->GetCollectionID( "IOTrackerCollection" );
    if(colIdIOTr>=0){
      IOTrHC=dynamic_cast<IOTrackerHitsCollection *>( HCE->GetHC( colIdIOTr ) );
      if(IOTrHC) nhIOTr=IOTrHC->entries();
    }
    G4int colIdSOTr=SDMan->GetCollectionID( "SOTrackerCollection" );
    if(colIdSOTr>=0){
      SOTrHC=dynamic_cast<SOTrackerHitsCollection *>( HCE->GetHC( colIdSOTr ) );
      if(SOTrHC) nhSOTr=SOTrHC->entries();
    }
    G4int colIdTof=SDMan->GetCollectionID( "TofCollection" ); 
    if(colIdTof>=0){
      TofHC=dynamic_cast<TofHitsCollection *>( HCE->GetHC( colIdTof ) ); 
      if(TofHC) nhTof=TofHC->entries();
    }
    G4int colIdITof=SDMan->GetCollectionID( "ITofCollection" ); 
    if(colIdITof>=0){
      ITofHC=dynamic_cast<ITofHitsCollection *>( HCE->GetHC( colIdITof ) ); 
      if(ITofHC) nhITof=ITofHC->entries();
    }
    G4int colIdPAD=SDMan->GetCollectionID( "PADCollection" );
    if(colIdPAD>=0){
      PADHC=dynamic_cast<PADHitsCollection *>( HCE->GetHC( colIdPAD ) );
      if(PADHC) nhPAD=PADHC->entries();
    }
    G4int colIdRICH=SDMan->GetCollectionID( "RICHCollection" );
    if(colIdRICH>=0){
      RICHHC=dynamic_cast<RICHHitsCollection *>( HCE->GetHC( colIdRICH ) );
      if(RICHHC) nhRICH=RICHHC->entries();
    }
    G4int colIdPID1=SDMan->GetCollectionID( "PID1Collection" );
    if(colIdPID1>=0){
      PID1HC=dynamic_cast<PID1HitsCollection *>( HCE->GetHC( colIdPID1 ) );
      if(PID1HC) nhPID1=PID1HC->entries();
    }
    G4int colIdPID2=SDMan->GetCollectionID( "PID2Collection" );
    if(colIdPID2>=0){
      PID2HC=dynamic_cast<PID2HitsCollection *>( HCE->GetHC( colIdPID2 ) );
      if(PID2HC) nhPID2=PID2HC->entries();
    }
    G4int colIdMF=SDMan->GetCollectionID( "MuonFilterCollection" );
    if(colIdMF>=0){
      MFHC=dynamic_cast<MuonFilterHitsCollection *>( HCE->GetHC( colIdMF ) );
      if(MFHC) nhMF=MFHC->entries();
    }
    G4int colIdVD=SDMan->GetCollectionID( "VDCollection" );
    if(colIdVD>=0){
      VDHC=dynamic_cast<VDHitsCollection *>( HCE->GetHC( colIdVD ) );
      if(VDHC) nhVD=VDHC->entries();
    }
  }

  std::ios::fmtflags oldFlags = ost.flags();
  std::size_t preSiz = ost.precision();
  ost.setf( std::ios::fixed );
  ost << std::setw(12) << anEvent->GetEventID()+1
      << std::endl;

  // Primary Information
  if( pInfo_ ){
    ost << PrimInfoF << std::endl; 
    //IP poimt
    ost.precision(5);
    ost << std::setw(12) << pInfo_->x/mm 
	<< std::setw(12) << pInfo_->y/mm 
	<< std::setw(12) << pInfo_->z/mm;
    //	<< std::endl;
    //Beam
    ost << std::setw(12) << pInfo_->pbeam/GeV
	<< std::setw(12) << pInfo_->ubeam 
	<< std::setw(12) << pInfo_->vbeam;
    //	<< std::endl;
    //Measured beam mom
    G4bool flag_VD1 = false;
    for( G4int i=0; i<nhVD; i++ ){
      VDHit *aHit = (*VDHC)[i];
      G4int layer = aHit->GetLayerID(); 
      G4double p = aHit->GetLMom().mag();
      if( layer==1 ){
    	ost << std::setw(12) << p/GeV; 
    	//	    << std::endl;
	flag_VD1=true;
      }
    }
    if( !flag_VD1 ) ost << std::setw(12) << 0.0; 
    //Reaction: Scattered (D*)
    ost << std::setw(12) << pInfo_->m1/GeV
	<< std::setw(12) << pInfo_->p1/GeV
	<< std::setw(12) << pInfo_->theta1/degree
	<< std::setw(12) << pInfo_->phi1/degree
	<< std::setw(12) << pInfo_->theta1cm/degree
	<< std::setw(12) << pInfo_->phi1cm/degree;
    //	<< std::endl;
    //Reaction: Produced (Yc*)
    ost << std::setw(12) << pInfo_->m2/GeV
	<< std::setw(12) << pInfo_->p2/GeV
	<< std::setw(12) << pInfo_->theta2/degree
	<< std::setw(12) << pInfo_->phi2/degree
	<< std::setw(12) << pInfo_->theta2cm/degree
	<< std::setw(12) << pInfo_->phi2cm/degree
	<< std::endl;
  }
  ost << PrimInfoFEnd << std::endl; 

  /////////////////////////////////////////////
  ////Full tracking output
  ////Tracks are not selected.
  /////////////////////////////////////////////

  ost << FullTrackF << std::endl; 
  ost << FullTrackTF << std::endl; 

  //Tracker Information
  //BFT
  for( G4int i=0; i<nhBTr; i++ ){
    BTrackerHit *aHit = (*BTrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  //SFT
  for( G4int i=0; i<nhSITr; i++ ){
    SITrackerHit *aHit = (*SITrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  //IT1
  for( G4int i=0; i<nhIITr; i++ ){
    IITrackerHit *aHit = (*IITrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  //IT2
  for( G4int i=0; i<nhIOTr; i++ ){
    IOTrackerHit *aHit = (*IOTrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  //AFT
  for( G4int i=0; i<nhAOTr; i++ ){
    AOTrackerHit *aHit = (*AOTrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  //ST1&2
  for( G4int i=0; i<nhSOTr; i++ ){
    SOTrackerHit *aHit = (*SOTrHC)[i];
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID() 
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::endl;
  }
  ost << FullTrackTFEnd << std::endl; 

  //Counter information
  ost << FullTrackCF << std::endl; 
  //T0
  for( G4int i=0; i<nhT0; i++ ){
    T0Hit *aHit = (*T0HC)[i];
    ost << geomMan.GetDetectorId("T0")
	<< std::setw(3) << nhT0 << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //Tof
  for( G4int i=0; i<nhTof; i++ ){
    TofHit *aHit = (*TofHC)[i];
    ost << geomMan.GetDetectorId("TOF")
	<< std::setw(3) << nhTof << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //ITof
  for( G4int i=0; i<nhITof; i++ ){
    ITofHit *aHit = (*ITofHC)[i];
    ost << geomMan.GetDetectorId("ITOFR")
	<< std::setw(3) << nhITof << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //PAD
  for( G4int i=0; i<nhPAD; i++ ){
    PADHit *aHit = (*PADHC)[i];
    ost << geomMan.GetDetectorId("PADU")
	<< std::setw(3) << nhPAD << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //RICH
  for( G4int i=0; i<nhRICH; i++ ){
    RICHHit *aHit = (*RICHHC)[i];
    ost << geomMan.GetDetectorId("RICH")
	<< std::setw(3) << nhRICH << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //PID1
  for( G4int i=0; i<nhPID1; i++ ){
    PID1Hit *aHit = (*PID1HC)[i];
    ost << geomMan.GetDetectorId("PID1R")
	<< std::setw(3) << nhPID1 << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //PID2
  for( G4int i=0; i<nhPID2; i++ ){
    PID2Hit *aHit = (*PID2HC)[i];
    ost << geomMan.GetDetectorId("PID2R")
	<< std::setw(3) << nhPID2 << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //MF
  for( G4int i=0; i<nhMF; i++ ){
    MuonFilterHit *aHit = (*MFHC)[i];
    ost << geomMan.GetDetectorId("MFU")
	<< std::setw(3) << nhMF << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  //VD
  for( G4int i=0; i<nhVD; i++ ){
    VDHit *aHit = (*VDHC)[i];
    ost << geomMan.GetDetectorId("SpecVp1")
	<< std::setw(3) << nhVD << std::endl;
    ost.precision(5);
    ost << std::setw(12) << aHit->GetLayerID()
  	<< std::setw(12) << aHit->GetSegmentID()
  	<< std::setw(12) << aHit->GetTime()/ns
  	<< std::setw(12) << aHit->GetEdep()/MeV
  	<< std::setw(12) << aHit->GetPathLength()/mm
  	<< std::setw(12) << aHit->GetLMom().mag()/GeV
  	<< std::setw(12) << aHit->GetXLocal()/mm
  	<< std::setw(12) << aHit->GetYLocal()/mm
  	<< std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    ost.precision(8);
    ost	<< std::setw(14) << aHit->GetBeta()
  	<< std::endl;
  }
  ost << FullTrackCFEnd << std::endl; 
  ost << FullTrackFEnd << std::endl; 


  /////////////////////////////////////////////
  ////Simple tracking output
  ////Good tracks are already selected.
  /////////////////////////////////////////////

  ost << SimpleTrackF << std::endl; 
  for( G4int i=0; i<nhTgt; i++ ){
    TargetHit *aHittgt=(*TgtHC)[i];
    G4int id=aHittgt->GetTrackNo();
    G4int Pid=PIDFunc(aHittgt->GetHitParticleName(0));
    G4ThreeVector mom=aHittgt->GetLMom();
    G4ThreeVector Mom(mom.x(), mom.y(), mom.z());
    G4double xl=(-1.)*aHittgt->GetXLocal(), yl=(-1.)*aHittgt->GetYLocal();
    G4int type=0;

    //T0 hit analysis
    for( G4int j=0; j<nhT0; j++ ){
      T0Hit *aHit=(*T0HC)[j];
      G4int trackNo=aHit->GetTrackNo();

      if( id==trackNo ){
  	//Beam information
  	ost << SimpleTrackBF << std::endl; 
  	ost.precision(1);
  	ost << id << std::endl;
	type=1;
  	ost << type << std::endl;

  	//Generated information
  	ost << SimpleTrackGF << std::endl;
  	ost.precision(5);
  	ost << std::setw(12) << Mom.mag()/GeV
  	    << std::setw(12) << (-1.)*Mom.x()/GeV
  	    << std::setw(12) << (-1.)*Mom.y()/GeV
  	    << std::setw(12) << (-1.)*Mom.z()/GeV
  	    << std::setw(12) << Pid
  	    << std::setw(12) << xl/mm
  	    << std::setw(12) << yl/mm
  	    << std::endl;
  	ost << SimpleTrackGFEnd << std::endl;

  	//Counter information
  	ost << SimpleTrackCF << std::endl;
	//T0
  	ost << geomMan.GetDetectorId("T0") << std::endl;
  	ost.precision(5);
  	ost << std::setw(12) << aHit->GetLayerID()
  	    << std::setw(12) << aHit->GetSegmentID()
  	    << std::setw(12) << aHit->GetTime()/ns
  	    << std::setw(12) << aHit->GetEdep()/MeV
  	    << std::setw(12) << aHit->GetPathLength()/mm
  	    << std::setw(12) << aHit->GetLMom().mag()/GeV
  	    << std::setw(12) << aHit->GetXLocal()/mm
  	    << std::setw(12) << aHit->GetYLocal()/mm
  	    << std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
  	ost.precision(8);
  	ost << std::setw(14) << aHit->GetBeta()
  	    << std::endl;
  	ost << SimpleTrackCFEnd << std::endl;	

	//Tracker information
  	ost << SimpleTrackTF << std::endl;	
  	//BFT
  	for( G4int k=0; k<nhBTr; k++ ){
  	  BTrackerHit *aHit1=(*BTrHC)[k];
  	  G4int trackNo1=aHit1->GetTrackNo();
  	  if( trackNo==trackNo1 ){
  	    ost.precision(5);
  	    ost << std::setw(12) << aHit1->GetLayerID() 
  		<< std::setw(12) << aHit1->GetXLocal()/mm
  		<< std::setw(12) << aHit1->GetYLocal()/mm
  		<< std::endl;
  	  }
  	}
  	ost << SimpleTrackTFEnd << std::endl;	
	ost << SimpleTrackBFEnd << std::endl; 
      }
    }

    //TOF hit analysis
    for( G4int j=0; j<nhTof; j++ ){
      TofHit *aHit=(*TofHC)[j];
      G4int trackNo=aHit->GetTrackNo();
	
      if( id==trackNo ){
    	//Scattered particle information
    	ost << SimpleTrackSF << std::endl; 
    	ost.precision(1);
    	ost << id << std::endl;
	type=2;
  	ost << type << std::endl;

    	//Generated information
    	ost << SimpleTrackGF << std::endl;
    	ost.precision(5);
    	ost << std::setw(12) << Mom.mag()/GeV
    	    << std::setw(12) << (-1.)*Mom.x()/GeV
    	    << std::setw(12) << (-1.)*Mom.y()/GeV
    	    << std::setw(12) << (-1.)*Mom.z()/GeV
    	    << std::setw(12) << Pid
    	    << std::setw(12) << xl/mm
    	    << std::setw(12) << yl/mm
    	    << std::endl;
  	ost << SimpleTrackGFEnd << std::endl;

    	//Counter information
	ost << SimpleTrackCF << std::endl;
	//TOF
    	ost << geomMan.GetDetectorId("TOF") << std::endl;
    	ost.precision(5);
    	ost << std::setw(12) << aHit->GetLayerID()
    	    << std::setw(12) << aHit->GetSegmentID()
    	    << std::setw(12) << aHit->GetTime()/ns
    	    << std::setw(12) << aHit->GetEdep()/MeV
    	    << std::setw(12) << aHit->GetPathLength()/mm
    	    << std::setw(12) << aHit->GetLMom().mag()/GeV
    	    << std::setw(12) << aHit->GetXLocal()/mm
    	    << std::setw(12) << aHit->GetYLocal()/mm
    	    << std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    	ost.precision(8);
    	ost << std::setw(14) << aHit->GetBeta()
    	    << std::endl;

    	//RICH
    	for( G4int k=0; k<nhRICH; k++ ){
    	  RICHHit *aHit1=(*RICHHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ){
    	    //RICH information
    	    ost << geomMan.GetDetectorId("RICH") << std::endl;
    	    ost.precision(5);
    	    ost << std::setw(12) << aHit1->GetLayerID()
    		<< std::setw(12) << aHit1->GetSegmentID()
    		<< std::setw(12) << aHit1->GetTime()/ns
    		<< std::setw(12) << aHit1->GetEdep()/MeV
    		<< std::setw(12) << aHit1->GetPathLength()/mm
    		<< std::setw(12) << aHit1->GetLMom().mag()/GeV
    		<< std::setw(12) << aHit1->GetXLocal()/mm
    		<< std::setw(12) << aHit1->GetYLocal()/mm
    		<< std::setw(12) << PIDFunc(aHit1->GetHitParticleName(0));
    	    ost.precision(8);
    	    ost << std::setw(14) << aHit1->GetBeta()
    		<< std::endl;
	    
    	  }
    	}

    	//PID1 ID
    	for( G4int k=0; k<nhPID1; k++ ){
    	  PID1Hit *aHit1=(*PID1HC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ){
    	    //PID1 information
    	    ost << geomMan.GetDetectorId("PID1R") << std::endl;
    	    ost.precision(5);
    	    ost << std::setw(12) << aHit1->GetLayerID()
    		<< std::setw(12) << aHit1->GetSegmentID()
    		<< std::setw(12) << aHit1->GetTime()/ns
    		<< std::setw(12) << aHit1->GetEdep()/MeV
    		<< std::setw(12) << aHit1->GetPathLength()/mm
    		<< std::setw(12) << aHit1->GetLMom().mag()/GeV
    		<< std::setw(12) << aHit1->GetXLocal()/mm
    		<< std::setw(12) << aHit1->GetYLocal()/mm
    		<< std::setw(12) << PIDFunc(aHit1->GetHitParticleName(0));
    	    ost.precision(8);
    	    ost << std::setw(14) << aHit1->GetBeta()
    		<< std::endl;
    	  }
	}
	ost << SimpleTrackCFEnd << std::endl;

	//Tracker information
	ost << SimpleTrackTF << std::endl;	
    	//SFT	
	for( G4int k=0; k<nhSITr; k++ ){
    	  SITrackerHit *aHit1=(*SITrHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ){
    	    //SFT information
    	    ost.precision(5);
    	    ost << std::setw(12) << aHit1->GetLayerID() 
    		<< std::setw(12) << aHit1->GetXLocal()/mm
    		<< std::setw(12) << aHit1->GetYLocal()/mm
    		<< std::endl;
    	  }
    	}

    	//IT1
    	for( G4int k=0; k<nhIITr; k++ ){
    	  IITrackerHit *aHit1=(*IITrHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ){
    	    //IT1 information
    	    ost.precision(5);
    	    ost << std::setw(12) << aHit1->GetLayerID() 
    		<< std::setw(12) << aHit1->GetXLocal()/mm
    		<< std::setw(12) << aHit1->GetYLocal()/mm
    		<< std::endl;
    	  }
    	}

    	//ST1&2
    	for( G4int k=0; k<nhSOTr; k++ ){
    	  SOTrackerHit *aHit1=(*SOTrHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ){
    	    //ST1&2 information
    	    ost.precision(5);
    	    ost << std::setw(12) << aHit1->GetLayerID() 
    		<< std::setw(12) << aHit1->GetXLocal()/mm
    		<< std::setw(12) << aHit1->GetYLocal()/mm
    		<< std::endl;
    	  }
    	}
	ost << SimpleTrackTFEnd << std::endl;	
	ost << SimpleTrackSFEnd << std::endl; 
      }
    }

    //ITOF hit analysis
    for( G4int j=0; j<nhITof; j++ ){
      ITofHit *aHit=(*ITofHC)[j];
      G4int trackNo=aHit->GetTrackNo();
	
      if( id==trackNo ){
    	//TOF ID
    	G4bool flag_tof =false;
    	for( G4int k=0; k<nhTof; k++ ){
    	  TofHit *aHit1=(*TofHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ) flag_tof =true;
    	}

    	if( !flag_tof ){
    	  //Scattered particle information
	  ost << SimpleTrackSF << std::endl; 
    	  ost.precision(1);
    	  ost << id << std::endl;
	  type=3;
	  ost << type << std::endl;

    	  //Generated information
    	  ost << SimpleTrackGF << std::endl;
    	  ost.precision(5);
    	  ost << std::setw(12) << Mom.mag()/GeV
    	      << std::setw(12) << (-1.)*Mom.x()/GeV
    	      << std::setw(12) << (-1.)*Mom.y()/GeV
    	      << std::setw(12) << (-1.)*Mom.z()/GeV
    	      << std::setw(12) << Pid
    	      << std::setw(12) << xl/mm
    	      << std::setw(12) << yl/mm
    	      << std::endl;
    	  ost << SimpleTrackGFEnd << std::endl;

    	  //Counter information
	  ost << SimpleTrackCF << std::endl;
	  //ITOF
	  ost << geomMan.GetDetectorId("ITOFR") << std::endl;
    	  ost.precision(5);
    	  ost << std::setw(12) << aHit->GetLayerID()
    	      << std::setw(12) << aHit->GetSegmentID()
    	      << std::setw(12) << aHit->GetTime()/ns
    	      << std::setw(12) << aHit->GetEdep()/MeV
    	      << std::setw(12) << aHit->GetPathLength()/mm
    	      << std::setw(12) << aHit->GetLMom().mag()/GeV
    	      << std::setw(12) << aHit->GetXLocal()/mm
    	      << std::setw(12) << aHit->GetYLocal()/mm
    	      << std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    	  ost.precision(8);
    	  ost << std::setw(14) << aHit->GetBeta()
    	      << std::endl;

    	  //PID1 ID
    	  for( G4int k=0; k<nhPID1; k++ ){
    	    PID1Hit *aHit1=(*PID1HC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      //PID1 information
	      ost << geomMan.GetDetectorId("PID1R") << std::endl;
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID()
    		  << std::setw(12) << aHit1->GetSegmentID()
    		  << std::setw(12) << aHit1->GetTime()/ns
    		  << std::setw(12) << aHit1->GetEdep()/MeV
    		  << std::setw(12) << aHit1->GetPathLength()/mm
    		  << std::setw(12) << aHit1->GetLMom().mag()/GeV
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::setw(12) << PIDFunc(aHit1->GetHitParticleName(0));
    	      ost.precision(8);
    	      ost << std::setw(14) << aHit1->GetBeta()
    		  << std::endl;
    	    }
    	  }

    	  //PID2 ID
    	  for( G4int k=0; k<nhPID2; k++ ){
    	    PID2Hit *aHit1=(*PID2HC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      //PID2 information
	      ost << geomMan.GetDetectorId("PID2R") << std::endl;
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID()
    		  << std::setw(12) << aHit1->GetSegmentID()
    		  << std::setw(12) << aHit1->GetTime()/ns
    		  << std::setw(12) << aHit1->GetEdep()/MeV
    		  << std::setw(12) << aHit1->GetPathLength()/mm
    		  << std::setw(12) << aHit1->GetLMom().mag()/GeV
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::setw(12) << PIDFunc(aHit1->GetHitParticleName(0));
    	      ost.precision(8);
    	      ost << std::setw(14) << aHit1->GetBeta()
    		  << std::endl;
    	    }
    	  }
	  ost << SimpleTrackCFEnd << std::endl;
	
	  //Tracker infotmation
	  ost << SimpleTrackTF << std::endl;	
    	  //SFT
    	  for( G4int k=0; k<nhSITr; k++ ){
    	    SITrackerHit *aHit1=(*SITrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }

    	  //IT1
    	  for( G4int k=0; k<nhIITr; k++ ){
    	    IITrackerHit *aHit1=(*IITrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }

    	  //IT2
    	  for( G4int k=0; k<nhIOTr; k++ ){
    	    IOTrackerHit *aHit1=(*IOTrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }

    	  //ST1&2
    	  for( G4int k=0; k<nhSOTr; k++ ){
    	    SOTrackerHit *aHit1=(*SOTrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }

    	  //AFT
    	  for( G4int k=0; k<nhAOTr; k++ ){
    	    AOTrackerHit *aHit1=(*AOTrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      //AFT information
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }
	  ost << SimpleTrackTFEnd << std::endl;	
	  ost << SimpleTrackSFEnd << std::endl; 
    	}
      }
    }

    //PAD hit analysis
    for( G4int j=0; j<nhPAD; j++ ){
      PADHit *aHit=(*PADHC)[j];
      G4int trackNo=aHit->GetTrackNo();

      if( id==trackNo ){
    	//TOF ID
    	G4bool flag_tof =false;
    	for( G4int k=0; k<nhTof; k++ ){
    	  TofHit *aHit1=(*TofHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ) flag_tof =true;
    	}
    	//ITOF ID
    	G4bool flag_itof =false;
    	for( G4int k=0; k<nhITof; k++ ){
    	  ITofHit *aHit1=(*ITofHC)[k];
    	  G4int trackNo1=aHit1->GetTrackNo();
    	  if( trackNo==trackNo1 ) flag_itof =true;
    	}

    	if( !flag_tof && !flag_itof ){
    	  //Scattered particle information
	  ost << SimpleTrackSF << std::endl; 
    	  ost.precision(1);
    	  ost << id << std::endl;
	  type=4;
	  ost << type << std::endl;

    	  //Generated information
	  ost << SimpleTrackGF << std::endl;
    	  ost.precision(5);
    	  ost << std::setw(12) << Mom.mag()/GeV
    	      << std::setw(12) << (-1.)*Mom.x()/GeV
    	      << std::setw(12) << (-1.)*Mom.y()/GeV
    	      << std::setw(12) << (-1.)*Mom.z()/GeV
    	      << std::setw(12) << Pid
    	      << std::setw(12) << xl/mm
    	      << std::setw(12) << yl/mm
    	      << std::endl;
	  ost << SimpleTrackGFEnd << std::endl;

	  //Counter information
	  ost << SimpleTrackCF << std::endl;
    	  //PAD
	  ost << geomMan.GetDetectorId("PADU") << std::endl;
    	  ost.precision(5);
    	  ost << std::setw(12) << aHit->GetLayerID()
    	      << std::setw(12) << aHit->GetSegmentID()
    	      << std::setw(12) << aHit->GetTime()/ns
    	      << std::setw(12) << aHit->GetEdep()/MeV
    	      << std::setw(12) << aHit->GetPathLength()/mm
    	      << std::setw(12) << aHit->GetLMom().mag()/GeV
    	      << std::setw(12) << aHit->GetXLocal()/mm
    	      << std::setw(12) << aHit->GetYLocal()/mm
    	      << std::setw(12) << PIDFunc(aHit->GetHitParticleName(0));
    	  ost.precision(8);
    	  ost << std::setw(14) << aHit->GetBeta()
    	      << std::endl;
	  ost << SimpleTrackCFEnd << std::endl;

	  //Tracker infotmation
	  ost << SimpleTrackTF << std::endl;	
    	  //SFT
    	  for( G4int k=0; k<nhSITr; k++ ){
    	    SITrackerHit *aHit1=(*SITrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }

    	  //IT1
    	  for( G4int k=0; k<nhIITr; k++ ){
    	    IITrackerHit *aHit1=(*IITrHC)[k];
    	    G4int trackNo1=aHit1->GetTrackNo();
    	    if( trackNo==trackNo1 ){
    	      //IT1 information
    	      ost.precision(5);
    	      ost << std::setw(12) << aHit1->GetLayerID() 
    		  << std::setw(12) << aHit1->GetXLocal()/mm
    		  << std::setw(12) << aHit1->GetYLocal()/mm
    		  << std::endl;
    	    }
    	  }
	  ost << SimpleTrackTFEnd << std::endl;	
	  ost << SimpleTrackSFEnd << std::endl; 
    	}
      }
    }
  }
  ost << SimpleTrackFEnd << std::endl;

  ost << EventEnd << std::endl;
  ost.flags( oldFlags );
  ost.precision( preSiz );
}


//////Histograms
void AnalysisSpec::InitializeEvent( void )
{
  //Primary
  event.priposx = -9999.0; 
  event.priposy = -9999.0; 
  event.priposz = -9999.0;
  event.prim1 = -9999.0; 
  event.prim2 = -9999.0;
  event.prip1 = -9999.0; 
  event.prip2 = -9999.0;
  event.prip1cm = -9999.0; 
  event.prip2cm = -9999.0;
  event.pritheta1 = -9999.0; 
  event.pritheta2 = -9999.0;
  event.pritheta1cm = -9999.0; 
  event.pritheta2cm = -9999.0;
  event.priphi1 = -9999.0; 
  event.priphi2 = -9999.0;
  event.priphi1cm = -9999.0; 
  event.priphi2cm = -9999.0;
  event.pricos1 = -9999.0; 
  event.pricos2 = -9999.0;
  event.pricos1cm = -9999.0; 
  event.pricos2cm = -9999.0;
  event.priu1 = -9999.0; 
  event.priu2 = -9999.0;
  event.priu1cm = -9999.0; 
  event.priu2cm = -9999.0;
  event.priv1 = -9999.0; 
  event.priv2 = -9999.0;
  event.priv1cm = -9999.0; 
  event.priv2cm = -9999.0;
  event.pribeta1 = -9999.0; 
  event.pribeta2 = -9999.0;
  event.pbeam = -9999.0; 
  event.ubeam = -9999.0;
  event.vbeam = -9999.0;
  event.prip3 = -9999.0; 
  event.priu3 = -9999.0;
  event.priv3 = -9999.0;

  ConfMan *confMan = ConfMan::GetConfManager();
  G4bool counterHits = confMan->DetRootFill();

  if( counterHits ){
    //Target
    event.tgtnhits = -1;
    event.tgtid.clear();
    event.tgtpid.clear();
    event.tgtedep.clear();
    event.tgttime.clear();
    event.tgtp.clear();
    event.tgttheta.clear();
    event.tgtphi.clear();
    event.tgtposx.clear();
    event.tgtposy.clear();
    event.tgtpath.clear();
    event.tgtbeta.clear();
    event.tgtmass.clear();
    
    //T0
    event.t0nhits = -1;
    event.t0layer.clear();
    event.t0seg.clear();
    event.t0id.clear();
    event.t0pid.clear();
    event.t0edep.clear();
    event.t0time.clear();
    event.t0p.clear();
    event.t0theta.clear();
    event.t0phi.clear();
    event.t0posx.clear();
    event.t0posy.clear();
    event.t0path.clear();
    event.t0beta.clear();
    event.t0mass.clear();
    
    //BTracker  
    event.btrnhits = -1;
    event.btrlayer.clear();
    event.btrid.clear();
    event.btrpid.clear();
    event.btredep.clear();
    event.btrtime.clear();
    event.btrp.clear();
    event.btrtheta.clear();
    event.btrphi.clear();
    event.btrposx.clear();
    event.btrposy.clear();
    event.btrpath.clear();
    event.btrbeta.clear();
    event.btrmass.clear();

    //SITracker  
    event.sitrnhits = -1;
    event.sitrlayer.clear();
    event.sitrid.clear();
    event.sitrpid.clear();
    event.sitredep.clear();
    event.sitrtime.clear();
    event.sitrp.clear();
    event.sitrtheta.clear();
    event.sitrphi.clear();
    event.sitrposx.clear();
    event.sitrposy.clear();
    event.sitrpath.clear();
    event.sitrbeta.clear();
    event.sitrmass.clear();
    
    //AOTracker  
    event.aotrnhits = -1;
    event.aotrlayer.clear();
    event.aotrid.clear();
    event.aotrpid.clear();
    event.aotredep.clear();
    event.aotrtime.clear();
    event.aotrp.clear();
    event.aotrtheta.clear();
    event.aotrphi.clear();
    event.aotrposx.clear();
    event.aotrposy.clear();
    event.aotrpath.clear();
    event.aotrbeta.clear();
    event.aotrmass.clear();
    
    //IITracker  
    event.iitrnhits = -1;
    event.iitrlayer.clear();
    event.iitrid.clear();
    event.iitrpid.clear();
    event.iitredep.clear();
    event.iitrtime.clear();
    event.iitrp.clear();
    event.iitrtheta.clear();
    event.iitrphi.clear();
    event.iitrposx.clear();
    event.iitrposy.clear();
    event.iitrpath.clear();
    event.iitrbeta.clear();
    event.iitrmass.clear();
    
    //IOTracker  
    event.iotrnhits = -1;
    event.iotrlayer.clear();
    event.iotrid.clear();
    event.iotrpid.clear();
    event.iotredep.clear();
    event.iotrtime.clear();
    event.iotrp.clear();
    event.iotrtheta.clear();
    event.iotrphi.clear();
    event.iotrposx.clear();
    event.iotrposy.clear();
    event.iotrpath.clear();
    event.iotrbeta.clear();
    event.iotrmass.clear();
    
    //SOTracker  
    event.sotrnhits = -1;
    event.sotrlayer.clear();
    event.sotrid.clear();
    event.sotrpid.clear();
    event.sotredep.clear();
    event.sotrtime.clear();
    event.sotrp.clear();
    event.sotrtheta.clear();
    event.sotrphi.clear();
    event.sotrposx.clear();
    event.sotrposy.clear();
    event.sotrpath.clear();
    event.sotrbeta.clear();
    event.sotrmass.clear();
    
    //TOF
    event.tofnhits = -1;
    event.toflayer.clear();
    event.tofseg.clear();
    event.tofid.clear();
    event.tofpid.clear();
    event.tofedep.clear();
    event.toftime.clear();
    event.tofp.clear();
    event.toftheta.clear();
    event.tofphi.clear();
    event.tofposx.clear();
    event.tofposy.clear();
    event.tofpath.clear();
    event.tofbeta.clear();
    event.tofmass.clear();
    event.tofm2_0.clear();
    event.tofm2_1.clear();
    event.tofm2_2.clear();
    
    //ITOF
    event.itofnhits = -1;
    event.itoflayer.clear();
    event.itofseg.clear();
    event.itofid.clear();
    event.itofpid.clear();
    event.itofedep.clear();
    event.itoftime.clear();
    event.itofp.clear();
    event.itoftheta.clear();
    event.itofphi.clear();
    event.itofposx.clear();
    event.itofposy.clear();
    event.itofpath.clear();
    event.itofbeta.clear();
    event.itofmass.clear();
    event.itofm2_0.clear();
    event.itofm2_1.clear();
    event.itofm2_2.clear();

    //PAD
    event.padnhits = -1;
    event.padlayer.clear();
    event.padseg.clear();
    event.padid.clear();
    event.padpid.clear();
    event.padedep.clear();
    event.padtime.clear();
    event.padp.clear();
    event.padtheta.clear();
    event.padphi.clear();
    event.padposx.clear();
    event.padposy.clear();
    event.padpath.clear();
    event.padbeta.clear();
    event.padmass.clear();
    event.padm2_0.clear();
    event.padm2_1.clear();
    event.padm2_2.clear();
    
    //RICH
    event.richnhits = -1;
    event.richlayer.clear();
    event.richseg.clear();
    event.richid.clear();
    event.richpid.clear();
    event.richedep.clear();
    event.richtime.clear();
    event.richp.clear();
    event.richtheta.clear();
    event.richphi.clear();
    event.richposx.clear();
    event.richposy.clear();
    event.richpath.clear();
    event.richbeta.clear();
    event.richmass.clear();
    
    //PID1
    event.pid1nhits = -1;
    event.pid1layer.clear();
    event.pid1seg.clear();
    event.pid1id.clear();
    event.pid1pid.clear();
    event.pid1edep.clear();
    event.pid1time.clear();
    event.pid1p.clear();
    event.pid1theta.clear();
    event.pid1phi.clear();
    event.pid1posx.clear();
    event.pid1posy.clear();
    event.pid1path.clear();
    event.pid1beta.clear();
    event.pid1mass.clear();
    
    //PID2
    event.pid2nhits = -1;
    event.pid2layer.clear();
    event.pid2seg.clear();
    event.pid2id.clear();
    event.pid2pid.clear();
    event.pid2edep.clear();
    event.pid2time.clear();
    event.pid2p.clear();
    event.pid2theta.clear();
    event.pid2phi.clear();
    event.pid2posx.clear();
    event.pid2posy.clear();
    event.pid2path.clear();
    event.pid2beta.clear();
    event.pid2mass.clear();
    
    //MF
    event.mfnhits = -1;
    event.mflayer.clear();
    event.mfseg.clear();
    event.mfid.clear();
    event.mfpid.clear();
    event.mfedep.clear();
    event.mftime.clear();
    event.mfp.clear();
    event.mftheta.clear();
    event.mfphi.clear();
    event.mfposx.clear();
    event.mfposy.clear();
    event.mfpath.clear();
    event.mfbeta.clear();
    event.mfmass.clear();
    
    //VD
    event.vdnhits = -1;
    event.vdlayer.clear();
    event.vdid.clear();
    event.vdpid.clear();
    event.vdedep.clear();
    event.vdtime.clear();
    event.vdp.clear();
    event.vdtheta.clear();
    event.vdphi.clear();
    event.vdposx.clear();
    event.vdposy.clear();
    event.vdpath.clear();
    event.vdbeta.clear();
    event.vdmass.clear();
  }

  //PID
  event.ntKp = -1;
  event.ntKm = -1;
  event.ntPip = -1;
  event.ntPim = -1;
  event.ntP = -1;
  event.ntPb = -1;
  event.ntMup = -1;
  event.ntMum = -1;
  event.ntEp = -1;
  event.ntEm = -1;

  event.pkp.clear();
  event.ukp.clear();
  event.vkp.clear();
  event.thetakp.clear();
  event.phikp.clear();

  event.pkm.clear();
  event.ukm.clear();
  event.vkm.clear();
  event.thetakm.clear();
  event.phikm.clear();
  
  event.ppip.clear();
  event.upip.clear();
  event.vpip.clear();
  event.thetapip.clear();
  event.phipip.clear();

  event.ppim.clear();
  event.upim.clear();
  event.vpim.clear();
  event.thetapim.clear();
  event.phipim.clear();

  event.pp.clear();
  event.up.clear();
  event.vp.clear();
  event.thetap.clear();
  event.phip.clear();

  event.ppb.clear();
  event.upb.clear();
  event.vpb.clear();
  event.thetapb.clear();
  event.phipb.clear();

  event.pmup.clear();
  event.umup.clear();
  event.vmup.clear();
  event.thetamup.clear();
  event.phimup.clear();

  event.pmum.clear();
  event.umum.clear();
  event.vmum.clear();
  event.thetamum.clear();
  event.phimum.clear();
  
  event.pep.clear();
  event.uep.clear();
  event.vep.clear();
  event.thetaep.clear();
  event.phiep.clear();

  event.pem.clear();
  event.uem.clear();
  event.vem.clear();
  event.thetaem.clear();
  event.phiem.clear();

  //PID flag
  event.ntdet.clear();
  event.ntrich = -1;
  event.ntpid1 = -1;
  event.ntpid2 = -1;
  event.nttof = -1;
  event.ntitof = -1;
  event.ntpad = -1;
  event.ntmf = -1;

  event.fkp.clear();
  event.fkm.clear();
  event.fpip.clear();
  event.fpim.clear();
  event.fp.clear();
  event.fpb.clear();
  event.fmup.clear();
  event.fmum.clear();
  event.fep.clear();
  event.fem.clear();

  //Output
  for( G4int it=0; it<MaxHits; it++ ){
    event.pkpx[it]= -9999.0;
    event.pkpy[it]= -9999.0;
    event.pkpz[it]= -9999.0;

    event.pkmx[it]= -9999.0;
    event.pkmy[it]= -9999.0;
    event.pkmz[it]= -9999.0;

    event.ppipx[it]= -9999.0;
    event.ppipy[it]= -9999.0;
    event.ppipz[it]= -9999.0;

    event.ppimx[it]= -9999.0;
    event.ppimy[it]= -9999.0;
    event.ppimz[it]= -9999.0;

    event.ppx[it]= -9999.0;
    event.ppy[it]= -9999.0;
    event.ppz[it]= -9999.0;

    event.ppbx[it]= -9999.0;
    event.ppby[it]= -9999.0;
    event.ppbz[it]= -9999.0;

    event.pmupx[it]= -9999.0;
    event.pmupy[it]= -9999.0;
    event.pmupz[it]= -9999.0;

    event.pmumx[it]= -9999.0;
    event.pmumy[it]= -9999.0;
    event.pmumz[it]= -9999.0;

    event.pepx[it]= -9999.0;
    event.pepy[it]= -9999.0;
    event.pepz[it]= -9999.0;

    event.pemx[it]= -9999.0;
    event.pemy[it]= -9999.0;
    event.pemz[it]= -9999.0;
  }

  event.ntB = -1;
  for( int it=0; it<MaxHits; it++ ){
    event.pB[it]= -9999.0;
    event.uB[it]= -9999.0;
    event.vB[it]= -9999.0;
    event.pBx[it]= -9999.0;
    event.pBy[it]= -9999.0;
    event.pBz[it]= -9999.0;
  }

  for( int it=0; it<MaxHits; it++ ){
    event.flkp[it]= -1;
    event.flkm[it]= -1;
    event.flpip[it]= -1;
    event.flpim[it]= -1;
    event.flp[it]= -1;
    event.flpb[it]= -1;
    event.flmup[it]= -1;
    event.flmum[it]= -1;
    event.flep[it]= -1;
    event.flem[it]= -1;
  }

  return;
}

void AnalysisSpec::DefineHistograms( void )
{ 

  new TFile( filename_, "recreate" );
  fActive_=true;

  //Tree
  HBTree("tree","tree of Sks");
  TTree *tree = dynamic_cast<TTree *>(gFile->Get("tree"));
   
  //Primary
  tree->Branch("priposx", &event.priposx, "priposx/D"); 
  tree->Branch("priposy", &event.priposy, "priposy/D"); 
  tree->Branch("priposz", &event.priposz, "priposz/D"); 
  tree->Branch("prim1", &event.prim1, "prim1/D"); 
  tree->Branch("prim2", &event.prim2, "prim2/D"); 
  tree->Branch("prip1", &event.prip1, "prip1/D"); 
  tree->Branch("prip2", &event.prip2, "prip2/D"); 
  tree->Branch("prip1cm", &event.prip1cm, "prip1cm/D"); 
  tree->Branch("prip2cm", &event.prip2cm, "prip2cm/D"); 
  tree->Branch("pritheta1", &event.pritheta1, "pritheta1/D"); 
  tree->Branch("pritheta2", &event.pritheta2, "pritheta2/D"); 
  tree->Branch("pritheta1cm", &event.pritheta1cm, "pritheta1cm/D"); 
  tree->Branch("pritheta2cm", &event.pritheta2cm, "pritheta2cm/D"); 
  tree->Branch("priphi1", &event.priphi1, "priphi1/D"); 
  tree->Branch("priphi2", &event.priphi2, "priphi2/D"); 
  tree->Branch("priphi1cm", &event.priphi1cm, "priphi1cm/D"); 
  tree->Branch("priphi2cm", &event.priphi2cm, "priphi2cm/D"); 
  tree->Branch("pricos1", &event.pricos1, "pricos1/D"); 
  tree->Branch("pricos2", &event.pricos2, "pricos2/D"); 
  tree->Branch("pricos1cm", &event.pricos1cm, "pricos1cm/D"); 
  tree->Branch("pricos2cm", &event.pricos2cm, "pricos2cm/D"); 
  tree->Branch("priu1", &event.priu1, "priu1/D"); 
  tree->Branch("priu2", &event.priu2, "priu2/D"); 
  tree->Branch("priu1cm", &event.priu1cm, "priu1cm/D"); 
  tree->Branch("priu2cm", &event.priu2cm, "priu2cm/D"); 
  tree->Branch("priv1", &event.priv1, "priv1/D"); 
  tree->Branch("priv2", &event.priv2, "priv2/D"); 
  tree->Branch("priv1cm", &event.priv1cm, "priv1cm/D"); 
  tree->Branch("priv2cm", &event.priv2cm, "priv2cm/D"); 
  tree->Branch("pribeta1", &event.pribeta1, "pribeta1/D"); 
  tree->Branch("pribeta2", &event.pribeta2, "pribeta2/D"); 
  tree->Branch("pbeam", &event.pbeam, "pbeam/D"); 
  tree->Branch("ubeam", &event.ubeam, "ubeam/D"); 
  tree->Branch("vbeam", &event.vbeam, "vbeam/D"); 
  tree->Branch("prip3", &event.prip3, "prip3/D"); 
  tree->Branch("priu3", &event.priu3, "priu3/D"); 
  tree->Branch("priv3", &event.priv3, "priv3/D"); 

  ConfMan *confMan = ConfMan::GetConfManager();
  G4bool counterHits = confMan->DetRootFill();

  if( counterHits ){
    //Terget
    tree->Branch("tgtnhits", &event.tgtnhits);
    tree->Branch("tgtpid", &event.tgtpid);
    tree->Branch("tgtid", &event.tgtid);
    tree->Branch("tgtedep", &event.tgtedep);
    tree->Branch("tgttime", &event.tgttime);
    tree->Branch("tgtp", &event.tgtp);
    tree->Branch("tgttheta", &event.tgttheta);
    tree->Branch("tgtphi", &event.tgtphi);
    tree->Branch("tgtposx", &event.tgtposx);
    tree->Branch("tgtposy", &event.tgtposy);
    tree->Branch("tgtpath", &event.tgtpath);
    tree->Branch("tgtbeta", &event.tgtbeta);
    tree->Branch("tgtmass", &event.tgtmass);
    
    //T0
    tree->Branch("t0nhits", &event.t0nhits);
    tree->Branch("t0layer", &event.t0layer);
    tree->Branch("t0seg", &event.t0seg);
    tree->Branch("t0pid", &event.t0pid);
    tree->Branch("t0id", &event.t0id);
    tree->Branch("t0edep", &event.t0edep);
    tree->Branch("t0time", &event.t0time);
    tree->Branch("t0p", &event.t0p);
    tree->Branch("t0theta", &event.t0theta);
    tree->Branch("t0phi", &event.t0phi);
    tree->Branch("t0posx", &event.t0posx);
    tree->Branch("t0posy", &event.t0posy);
    tree->Branch("t0path", &event.t0path);
    tree->Branch("t0beta", &event.t0beta);
    tree->Branch("t0mass", &event.t0mass);
    
    //BTracker
    tree->Branch("btrnhits", &event.btrnhits);
    tree->Branch("btrlayer", &event.btrlayer);
    tree->Branch("btrpid", &event.btrpid);
    tree->Branch("btrid", &event.btrid);
    tree->Branch("btredep", &event.btredep);
    tree->Branch("btrtime", &event.btrtime);
    tree->Branch("btrp", &event.btrp);
    tree->Branch("btrtheta", &event.btrtheta);
    tree->Branch("btrphi", &event.btrphi);
    tree->Branch("btrposx", &event.btrposx);
    tree->Branch("btrposy", &event.btrposy);
    tree->Branch("btrpath", &event.btrpath);
    tree->Branch("btrbeta", &event.btrbeta);
    tree->Branch("btrmass", &event.btrmass);
    
    //SITracker
    tree->Branch("sitrnhits", &event.sitrnhits);
    tree->Branch("sitrlayer", &event.sitrlayer);
    tree->Branch("sitrpid", &event.sitrpid);
    tree->Branch("sitrid", &event.sitrid);
    tree->Branch("sitredep", &event.sitredep);
    tree->Branch("sitrtime", &event.sitrtime);
    tree->Branch("sitrp", &event.sitrp);
    tree->Branch("sitrtheta", &event.sitrtheta);
    tree->Branch("sitrphi", &event.sitrphi);
    tree->Branch("sitrposx", &event.sitrposx);
    tree->Branch("sitrposy", &event.sitrposy);
    tree->Branch("sitrpath", &event.sitrpath);
    tree->Branch("sitrbeta", &event.sitrbeta);
    tree->Branch("sitrmass", &event.sitrmass);

    //AOTracker
    tree->Branch("aotrnhits", &event.aotrnhits);
    tree->Branch("aotrlayer", &event.aotrlayer);
    tree->Branch("aotrpid", &event.aotrpid);
    tree->Branch("aotrid", &event.aotrid);
    tree->Branch("aotredep", &event.aotredep);
    tree->Branch("aotrtime", &event.aotrtime);
    tree->Branch("aotrp", &event.aotrp);
    tree->Branch("aotrtheta", &event.aotrtheta);
    tree->Branch("aotrphi", &event.aotrphi);
    tree->Branch("aotrposx", &event.aotrposx);
    tree->Branch("aotrposy", &event.aotrposy);
    tree->Branch("aotrpath", &event.aotrpath);
    tree->Branch("aotrbeta", &event.aotrbeta);
    tree->Branch("aotrmass", &event.aotrmass);

    //IITracker
    tree->Branch("iitrnhits", &event.iitrnhits);
    tree->Branch("iitrlayer", &event.iitrlayer);
    tree->Branch("iitrpid", &event.iitrpid);
    tree->Branch("iitrid", &event.iitrid);
    tree->Branch("iitredep", &event.iitredep);
    tree->Branch("iitrtime", &event.iitrtime);
    tree->Branch("iitrp", &event.iitrp);
    tree->Branch("iitrtheta", &event.iitrtheta);
    tree->Branch("iitrphi", &event.iitrphi);
    tree->Branch("iitrposx", &event.iitrposx);
    tree->Branch("iitrposy", &event.iitrposy);
    tree->Branch("iitrpath", &event.iitrpath);
    tree->Branch("iitrbeta", &event.iitrbeta);
    tree->Branch("iitrmass", &event.iitrmass);
    
    //IOTracker
    tree->Branch("iotrnhits", &event.iotrnhits);
    tree->Branch("iotrlayer", &event.iotrlayer);
    tree->Branch("iotrpid", &event.iotrpid);
    tree->Branch("iotrid", &event.iotrid);
    tree->Branch("iotredep", &event.iotredep);
    tree->Branch("iotrtime", &event.iotrtime);
    tree->Branch("iotrp", &event.iotrp);
    tree->Branch("iotrtheta", &event.iotrtheta);
    tree->Branch("iotrphi", &event.iotrphi);
    tree->Branch("iotrposx", &event.iotrposx);
    tree->Branch("iotrposy", &event.iotrposy);
    tree->Branch("iotrpath", &event.iotrpath);
    tree->Branch("iotrbeta", &event.iotrbeta);
    tree->Branch("iotrmass", &event.iotrmass);

    //SOTracker
    tree->Branch("sotrnhits", &event.sotrnhits);
    tree->Branch("sotrlayer", &event.sotrlayer);
    tree->Branch("sotrpid", &event.sotrpid);
    tree->Branch("sotrid", &event.sotrid);
    tree->Branch("sotredep", &event.sotredep);
    tree->Branch("sotrtime", &event.sotrtime);
    tree->Branch("sotrp", &event.sotrp);
    tree->Branch("sotrtheta", &event.sotrtheta);
    tree->Branch("sotrphi", &event.sotrphi);
    tree->Branch("sotrposx", &event.sotrposx);
    tree->Branch("sotrposy", &event.sotrposy);
    tree->Branch("sotrpath", &event.sotrpath);
    tree->Branch("sotrbeta", &event.sotrbeta);
    tree->Branch("sotrmass", &event.sotrmass);
    
    //TOF
    tree->Branch("tofnhits", &event.tofnhits);
    tree->Branch("toflayer", &event.toflayer);
    tree->Branch("tofseg", &event.tofseg);
    tree->Branch("tofpid", &event.tofpid);
    tree->Branch("tofid", &event.tofid);
    tree->Branch("tofedep", &event.tofedep);
    tree->Branch("toftime", &event.toftime);
    tree->Branch("tofp", &event.tofp);
    tree->Branch("toftheta", &event.toftheta);
    tree->Branch("tofphi", &event.tofphi);
    tree->Branch("tofposx", &event.tofposx);
    tree->Branch("tofposy", &event.tofposy);
    tree->Branch("tofpath", &event.tofpath);
    tree->Branch("tofbeta", &event.tofbeta);
    tree->Branch("tofmass", &event.tofmass);
    tree->Branch("tofm2_0", &event.tofm2_0);
    tree->Branch("tofm2_1", &event.tofm2_1);
    tree->Branch("tofm2_2", &event.tofm2_2);
    
    //ITOF
    tree->Branch("itofnhits", &event.itofnhits);
    tree->Branch("itoflayer", &event.itoflayer);
    tree->Branch("itofseg", &event.itofseg);
    tree->Branch("itofpid", &event.itofpid);
    tree->Branch("itofid", &event.itofid);
    tree->Branch("itofedep", &event.itofedep);
    tree->Branch("itoftime", &event.itoftime);
    tree->Branch("itofp", &event.itofp);
    tree->Branch("itoftheta", &event.itoftheta);
    tree->Branch("itofphi", &event.itofphi);
    tree->Branch("itofposx", &event.itofposx);
    tree->Branch("itofposy", &event.itofposy);
    tree->Branch("itofpath", &event.itofpath);
    tree->Branch("itofbeta", &event.itofbeta);
    tree->Branch("itofmass", &event.itofmass);
    tree->Branch("itofm2_0", &event.itofm2_0);
    tree->Branch("itofm2_1", &event.itofm2_1);
    tree->Branch("itofm2_2", &event.itofm2_2);

    //PAD
    tree->Branch("padnhits", &event.padnhits);
    tree->Branch("padlayer", &event.padlayer);
    tree->Branch("padseg", &event.padseg);
    tree->Branch("padpid", &event.padpid);
    tree->Branch("padid", &event.padid);
    tree->Branch("padedep", &event.padedep);
    tree->Branch("padtime", &event.padtime);
    tree->Branch("padp", &event.padp);
    tree->Branch("padtheta", &event.padtheta);
    tree->Branch("padphi", &event.padphi);
    tree->Branch("padposx", &event.padposx);
    tree->Branch("padposy", &event.padposy);
    tree->Branch("padpath", &event.padpath);
    tree->Branch("padbeta", &event.padbeta);
    tree->Branch("padmass", &event.padmass);
    tree->Branch("padm2_0", &event.padm2_0);
    tree->Branch("padm2_1", &event.padm2_1);
    tree->Branch("padm2_2", &event.padm2_2);
    
    //RICH
    tree->Branch("richnhits", &event.richnhits);
    tree->Branch("richlayer", &event.richlayer);
    tree->Branch("richseg", &event.richseg);
    tree->Branch("richpid", &event.richpid);
    tree->Branch("richid", &event.richid);
    tree->Branch("richedep", &event.richedep);
    tree->Branch("richtime", &event.richtime);
    tree->Branch("richp", &event.richp);
    tree->Branch("richtheta", &event.richtheta);
    tree->Branch("richphi", &event.richphi);
    tree->Branch("richposx", &event.richposx);
    tree->Branch("richposy", &event.richposy);
    tree->Branch("richpath", &event.richpath);
    tree->Branch("richbeta", &event.richbeta);
    tree->Branch("richmass", &event.richmass);
    
    //PID1
    tree->Branch("pid1nhits", &event.pid1nhits);
    tree->Branch("pid1layer", &event.pid1layer);
    tree->Branch("pid1seg", &event.pid1seg);
    tree->Branch("pid1pid", &event.pid1pid);
    tree->Branch("pid1id", &event.pid1id);
    tree->Branch("pid1edep", &event.pid1edep);
    tree->Branch("pid1time", &event.pid1time);
    tree->Branch("pid1p", &event.pid1p);
    tree->Branch("pid1theta", &event.pid1theta);
    tree->Branch("pid1phi", &event.pid1phi);
    tree->Branch("pid1posx", &event.pid1posx);
    tree->Branch("pid1posy", &event.pid1posy);
    tree->Branch("pid1path", &event.pid1path);
    tree->Branch("pid1beta", &event.pid1beta);
    tree->Branch("pid1mass", &event.pid1mass);
    
    //PID2
    tree->Branch("pid2nhits", &event.pid2nhits);
    tree->Branch("pid2layer", &event.pid2layer);
    tree->Branch("pid2seg", &event.pid2seg);
    tree->Branch("pid2pid", &event.pid2pid);
    tree->Branch("pid2id", &event.pid2id);
    tree->Branch("pid2edep", &event.pid2edep);
    tree->Branch("pid2time", &event.pid2time);
    tree->Branch("pid2p", &event.pid2p);
    tree->Branch("pid2theta", &event.pid2theta);
    tree->Branch("pid2phi", &event.pid2phi);
    tree->Branch("pid2posx", &event.pid2posx);
    tree->Branch("pid2posy", &event.pid2posy);
    tree->Branch("pid2path", &event.pid2path);
    tree->Branch("pid2beta", &event.pid2beta);
    tree->Branch("pid2mass", &event.pid2mass);
    
    //MF
    tree->Branch("mfnhits", &event.mfnhits);
    tree->Branch("mflayer", &event.mflayer);
    tree->Branch("mfseg", &event.mfseg);
    tree->Branch("mfpid", &event.mfpid);
    tree->Branch("mfid", &event.mfid);
    tree->Branch("mfedep", &event.mfedep);
    tree->Branch("mftime", &event.mftime);
    tree->Branch("mfp", &event.mfp);
    tree->Branch("mftheta", &event.mftheta);
    tree->Branch("mfphi", &event.mfphi);
    tree->Branch("mfposx", &event.mfposx);
    tree->Branch("mfposy", &event.mfposy);
    tree->Branch("mfpath", &event.mfpath);
    tree->Branch("mfbeta", &event.mfbeta);
    tree->Branch("mfmass", &event.mfmass);
    
    //VD
    tree->Branch("vdnhits", &event.vdnhits);
    tree->Branch("vdlayer", &event.vdlayer);
    tree->Branch("vdpid", &event.vdpid);
    tree->Branch("vdid", &event.vdid);
    tree->Branch("vdedep", &event.vdedep);
    tree->Branch("vdtime", &event.vdtime);
    tree->Branch("vdp", &event.vdp);
    tree->Branch("vdtheta", &event.vdtheta);
    tree->Branch("vdphi", &event.vdphi);
    tree->Branch("vdposx", &event.vdposx);
    tree->Branch("vdposy", &event.vdposy);
    tree->Branch("vdpath", &event.vdpath);
    tree->Branch("vdbeta", &event.vdbeta);
    tree->Branch("vdmass", &event.vdmass);
  }

  //PID
  tree->Branch("ntKp", &event.ntKp, "ntKp/I");
  tree->Branch("ntKm", &event.ntKm, "ntKm/I");
  tree->Branch("ntPip", &event.ntPip, "ntPip/I");
  tree->Branch("ntPim", &event.ntPim, "ntPim/I");
  tree->Branch("ntP", &event.ntP, "ntP/I");
  tree->Branch("ntPb", &event.ntPb, "ntPb/I");
  tree->Branch("ntMup", &event.ntMup, "ntMup/I");
  tree->Branch("ntMum", &event.ntMum, "ntMum/I");
  tree->Branch("ntEp", &event.ntEp, "ntEp/I");
  tree->Branch("ntEm", &event.ntEm, "ntEm/I");

  tree->Branch("pkp", &event.pkp);
  tree->Branch("ukp", &event.ukp);
  tree->Branch("vkp", &event.vkp);
  tree->Branch("thetakp", &event.thetakp);
  tree->Branch("phikp", &event.phikp);

  tree->Branch("pkm", &event.pkm);
  tree->Branch("ukm", &event.ukm);
  tree->Branch("vkm", &event.vkm);
  tree->Branch("thetakm", &event.thetakm);
  tree->Branch("phikm", &event.phikm);

  tree->Branch("ppip", &event.ppip);
  tree->Branch("upip", &event.upip);
  tree->Branch("vpip", &event.vpip);
  tree->Branch("thetapip", &event.thetapip);
  tree->Branch("phipip", &event.phipip);

  tree->Branch("ppim", &event.ppim);
  tree->Branch("upim", &event.upim);
  tree->Branch("vpim", &event.vpim);
  tree->Branch("thetapim", &event.thetapim);
  tree->Branch("phipim", &event.phipim);

  tree->Branch("pp", &event.pp);
  tree->Branch("up", &event.up);
  tree->Branch("vp", &event.vp);
  tree->Branch("thetap", &event.thetap);
  tree->Branch("phip", &event.phip);

  tree->Branch("ppb", &event.ppb);
  tree->Branch("upb", &event.upb);
  tree->Branch("vpb", &event.vpb);
  tree->Branch("thetapb", &event.thetapb);
  tree->Branch("phipb", &event.phipb);

  tree->Branch("pmup", &event.pmup);
  tree->Branch("umup", &event.umup);
  tree->Branch("vmup", &event.vmup);
  tree->Branch("thetamup", &event.thetamup);
  tree->Branch("phimup", &event.phimup);

  tree->Branch("pmum", &event.pmum);
  tree->Branch("umum", &event.umum);
  tree->Branch("vmum", &event.vmum);
  tree->Branch("thetamum", &event.thetamum);
  tree->Branch("phimum", &event.phimum);

  tree->Branch("pep", &event.pep);
  tree->Branch("uep", &event.uep);
  tree->Branch("vep", &event.vep);
  tree->Branch("thetaep", &event.thetaep);
  tree->Branch("phiep", &event.phiep);

  tree->Branch("pem", &event.pem);
  tree->Branch("uem", &event.uem);
  tree->Branch("vem", &event.vem);
  tree->Branch("thetaem", &event.thetaem);
  tree->Branch("phiem", &event.phiem);

  //PID flag
  tree->Branch("ntdet", &event.ntdet);
  tree->Branch("ntrich", &event.ntrich, "ntrich/I");
  tree->Branch("ntpid1", &event.ntpid1, "ntpid1/I");
  tree->Branch("ntpid2", &event.ntpid2, "ntpid2/I");
  tree->Branch("nttof", &event.nttof, "nttof/I");
  tree->Branch("ntitof", &event.ntitof, "ntitof/I");
  tree->Branch("ntpad", &event.ntpad, "ntpad/I");
  tree->Branch("ntmf", &event.ntmf, "ntmf/I");

  tree->Branch("fkp", &event.fkp);
  tree->Branch("fkm", &event.fkm);
  tree->Branch("fpip", &event.fpip);
  tree->Branch("fpim", &event.fpim);
  tree->Branch("fp", &event.fp);
  tree->Branch("fpb", &event.fpb);
  tree->Branch("fmup", &event.fmup);
  tree->Branch("fmum", &event.fmum);
  tree->Branch("fep", &event.fep);
  tree->Branch("fem", &event.fem);

  //Output
  tree->Branch("pkpx", event.pkpx, "pkpx[ntKp]/D");
  tree->Branch("pkpy", event.pkpy, "pkpy[ntKp]/D");
  tree->Branch("pkpz", event.pkpz, "pkpz[ntKp]/D");

  tree->Branch("pkmx", event.pkmx, "pkmx[ntKm]/D");
  tree->Branch("pkmy", event.pkmy, "pkmy[ntKm]/D");
  tree->Branch("pkmz", event.pkmz, "pkmz[ntKm]/D");

  tree->Branch("ppipx", event.ppipx, "ppipx[ntPip]/D");
  tree->Branch("ppipy", event.ppipy, "ppipy[ntPip]/D");
  tree->Branch("ppipz", event.ppipz, "ppipz[ntPip]/D");

  tree->Branch("ppimx", event.ppimx, "ppimx[ntPim]/D");
  tree->Branch("ppimy", event.ppimy, "ppimy[ntPim]/D");
  tree->Branch("ppimz", event.ppimz, "ppimz[ntPim]/D");

  tree->Branch("ppx", event.ppx, "ppx[ntP]/D");
  tree->Branch("ppy", event.ppy, "ppy[ntP]/D");
  tree->Branch("ppz", event.ppz, "ppz[ntP]/D");

  tree->Branch("ppbx", event.ppbx, "ppbx[ntPb]/D");
  tree->Branch("ppby", event.ppby, "ppby[ntPb]/D");
  tree->Branch("ppbz", event.ppbz, "ppbz[ntPb]/D");

  tree->Branch("pmupx", event.pmupx, "pmupx[ntMup]/D");
  tree->Branch("pmupy", event.pmupy, "pmupy[ntMup]/D");
  tree->Branch("pmupz", event.pmupz, "pmupz[ntMup]/D");

  tree->Branch("pmumx", event.pmumx, "pmumx[ntMum]/D");
  tree->Branch("pmumy", event.pmumy, "pmumy[ntMum]/D");
  tree->Branch("pmumz", event.pmumz, "pmumz[ntMum]/D");

  tree->Branch("pepx", event.pepx, "pepx[ntEp]/D");
  tree->Branch("pepy", event.pepy, "pepy[ntEp]/D");
  tree->Branch("pepz", event.pepz, "pepz[ntEp]/D");

  tree->Branch("pemx", event.pemx, "pemx[ntEm]/D");
  tree->Branch("pemy", event.pemy, "pemy[ntEm]/D");
  tree->Branch("pemz", event.pemz, "pemz[ntEm]/D");

  tree->Branch("ntB", &event.ntB, "ntB/I");
  tree->Branch("pB", event.pB, "pB[ntB]/D");
  tree->Branch("uB", event.uB, "uB[ntB]/D");
  tree->Branch("vB", event.vB, "vB[ntB]/D");
  tree->Branch("pBx", event.pBx, "pBx[ntB]/D");
  tree->Branch("pBy", event.pBy, "pBy[ntB]/D");
  tree->Branch("pBz", event.pBz, "pBz[ntB]/D");

  tree->Branch("flkp", event.flkp, "flkp[ntKp]/I");
  tree->Branch("flkm", event.flkm, "flkm[ntKm]/I");
  tree->Branch("flpip", event.flpip, "flpip[ntPip]/I");
  tree->Branch("flpim", event.flpim, "flpim[ntPim]/I");
  tree->Branch("flp", event.flp, "flp[ntP]/I");
  tree->Branch("flpb", event.flpb, "flp[ntPb]/I");
  tree->Branch("flmup", event.flmup, "flmup[ntMup]/I");
  tree->Branch("flmum", event.flmum, "flmum[ntMum]/I");
  tree->Branch("flep", event.flep, "flep[ntEp]/I");
  tree->Branch("flem", event.flem, "flem[ntEm]/I");
}
