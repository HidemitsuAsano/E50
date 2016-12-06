/*
  AnalysisSpec.cc

  K.Shirotori
*/

#include "AnalysisSpec.hh"
#include "PrimaryInfo.hh"
#include "DataType.hh"

#include "T0Hit.hh"
#include "SFTHit.hh"

#include "PidUtil.hh"
#include "ConfMan.hh"
#include "DCGeomMan.hh"

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

struct Event{
  //Primary info
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
 // G4int    pnpart;

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

  //SFT
  G4int    sftnhits;
  std::vector<G4int>    sftlayer, sftseg;
  std::vector<G4int>    sftpid, sftid;
  std::vector<G4double> sftedep, sfttime;
  std::vector<G4double> sftp;
  std::vector<G4double> sfttheta, sftphi;
  std::vector<G4double> sftposlx, sftposly;
  std::vector<G4double> sftposgx, sftposgy, sftposgz;
  std::vector<G4double> sftpath, sftbeta;
  std::vector<G4double> sftmass;
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
    G4double x=pInfo_->x/mm, y=pInfo_->y/mm, z=pInfo_->z/mm; // primary vertex info
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

  G4int nhT0=0, nhSFT=0;
  T0HitsCollection         *T0HC=0;
  SFTHitsCollection  *SFTHC=0;

  if(HCE){
    G4int colIdT0=SDMan->GetCollectionID( "T0Collection" );
    if(colIdT0>=0){
      T0HC=dynamic_cast<T0HitsCollection *>( HCE->GetHC( colIdT0 ) );
      if(T0HC) nhT0=T0HC->entries();
    }
    G4int colIdSFT=SDMan->GetCollectionID( "SFTCollection" );
    if(colIdSFT>=0){
      SFTHC=dynamic_cast<SFTHitsCollection *>( HCE->GetHC( colIdSFT ) );
      if(SFTHC) nhSFT=SFTHC->entries();
    }
  }

  ConfMan *confMan = ConfMan::GetConfManager();
  G4bool counterHits = confMan->DetRootFill();

  if( counterHits ){
    ///////////////////////////////////////////////////////////////////
    ////////Counter Information
    ///////////////////////////

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
      G4ThreeVector Mom(mom.x(), mom.y(), mom.z());
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
    
    //SFT
    event.sftnhits = nhSFT;
    for( G4int i=0; i<nhSFT; i++ ){
      SFTHit *aHit=(*SFTHC)[i];
      
      G4int layer=aHit->GetLayerID();
      G4int seg=aHit->GetSegmentID();
      G4int id=aHit->GetTrackNo();
      G4int Pid=PIDFunc(aHit->GetHitParticleName(0));
      G4double edep=aHit->GetEdep();
      G4double time=aHit->GetTime();
      G4ThreeVector mom=aHit->GetLMom();
      G4ThreeVector Mom(mom.x(), mom.y(), mom.z());
      G4double xl=aHit->GetXLocal(), yl=aHit->GetYLocal();
      G4double xg=aHit->GetXGlobal();
      G4double yg=aHit->GetYGlobal();
      G4double zg=aHit->GetZGlobal();
      G4double p=Mom.mag();    
      G4double theta=Mom.theta()*Rad2Deg;
      G4double phi=Mom.phi()*Rad2Deg;
      G4double path=aHit->GetPathLength();
      G4double beta=aHit->GetBeta();
      G4double mass=aHit->GetMass();
    
      //         G4cout<< i << " id=" << id << " Layer=" << layer
      //         	  << " Mass=" << mass 
      //     	  << " PID=" << Pid << ":" << aHit->GetHitParticleName(0) << G4endl;
      
      event.sftlayer.push_back(layer);
      event.sftseg.push_back(seg);
      event.sftpid.push_back(Pid);
      event.sftid.push_back(id);
      event.sftedep.push_back(edep/MeV);
      event.sfttime.push_back(time/ns);
      event.sftp.push_back(p/GeV);
      event.sfttheta.push_back(theta);
      event.sftphi.push_back(phi);
      event.sftposlx.push_back(xl/mm);
      event.sftposly.push_back(yl/mm);
      event.sftposgx.push_back(xg/mm);
      event.sftposgy.push_back(yg/mm);
      event.sftposgz.push_back(zg/mm);
      event.sftpath.push_back(path/mm);
      event.sftbeta.push_back(beta);//
      event.sftmass.push_back(mass);// particle mass
    }
  }    

  G4bool trigCond1=false, trigCond2=false;
  if( nhT0>=1 ){ trigCond1 = true;}
  
  trigCond1 = true;
  trigCond2 = true;

  if( DataFile_.is_open() && trigCond1 && trigCond2 )
    PrintHitsInformation( anEvent, DataFile_ );
  fTriggered = ( trigCond1 && trigCond2 );


  tree->Fill();  
  
  //////////////////////End
}

void AnalysisSpec::PrintHitsInformation( const G4Event *anEvent, std::ostream &ost ) const
{
  G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  ConfMan *confMan = ConfMan::GetConfManager();

  G4int nhT0=0, nhSFT=0;
  T0HitsCollection         *T0HC=0;
  SFTHitsCollection  *SFTHC=0;

  if(HCE){
    G4int colIdT0=SDMan->GetCollectionID( "T0Collection" );
    if(colIdT0>=0){
      T0HC=dynamic_cast<T0HitsCollection *>( HCE->GetHC( colIdT0 ) );
      if(T0HC) nhT0=T0HC->entries();
    }
    G4int colIdSFT=SDMan->GetCollectionID( "SFTCollection" );
    if(colIdSFT>=0){
      SFTHC=dynamic_cast<SFTHitsCollection *>( HCE->GetHC( colIdSFT ) );
      if(SFTHC) nhSFT=SFTHC->entries();
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
    // for( G4int i=0; i<nhVD; i++ ){
    //   VDHit *aHit = (*VDHC)[i];
    //   G4int layer = aHit->GetLayerID(); 
    //   G4double p = aHit->GetLMom().mag();
    //   if( layer==1 ){
    // 	ost << std::setw(12) << p/GeV; 
    // 	//	    << std::endl;
    // 	flag_VD1=true;
    //   }
    // }
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
  //SFT
  for( G4int i=0; i<nhSFT; i++ ){
    SFTHit *aHit = (*SFTHC)[i];

    //Type 0, simple detector
    if( aHit->GetEdep() > confMan->GetTrEdep() &&
	confMan->DetectorType()==0 ){
      ost.precision(5);
      ost << std::setw(12) << aHit->GetLayerID()
	  << std::setw(12) << aHit->GetXLocal()/mm
	  << std::setw(12) << aHit->GetYLocal()/mm
	  << std::endl;
    }
    //realistic detector
    if( aHit->GetEdep() > confMan->GetTrEdep() &&
	confMan->DetectorType()>=1 ){
      ost.precision(5);
      ost << std::setw(12) << aHit->GetLayerID() 
	  << std::setw(12) << aHit->GetSegmentID()
	  << std::endl;
    }
  }
  ost << FullTrackTFEnd << std::endl; 

  //Counter information
  ost << FullTrackCF << std::endl; 
  //T0
  for( G4int i=0; i<nhT0; i++ ){
    T0Hit *aHit = (*T0HC)[i];
    if( aHit->GetEdep() > confMan->GetCoEdep() ){
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
  }
  ost << FullTrackCFEnd << std::endl; 
  ost << FullTrackFEnd << std::endl; 

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

    //SFT  
    event.sftnhits = -1;
    event.sftlayer.clear();
    event.sftseg.clear();
    event.sftid.clear();
    event.sftpid.clear();
    event.sftedep.clear();
    event.sfttime.clear();
    event.sftp.clear();
    event.sfttheta.clear();
    event.sftphi.clear();
    event.sftposlx.clear();
    event.sftposly.clear();
    event.sftposgx.clear();
    event.sftposgy.clear();
    event.sftposgz.clear();
    event.sftpath.clear();
    event.sftbeta.clear();
    event.sftmass.clear();
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
    
    //SFT
    tree->Branch("sftnhits", &event.sftnhits);
    tree->Branch("sftlayer", &event.sftlayer);
    tree->Branch("sftseg", &event.sftseg);
    tree->Branch("sftpid", &event.sftpid);
    tree->Branch("sftid", &event.sftid);
    tree->Branch("sftedep", &event.sftedep);
    tree->Branch("sfttime", &event.sfttime);
    tree->Branch("sftp", &event.sftp);
    tree->Branch("sfttheta", &event.sfttheta);
    tree->Branch("sftphi", &event.sftphi);
    tree->Branch("sftposlx", &event.sftposlx);
    tree->Branch("sftposly", &event.sftposly);
    tree->Branch("sftposgx", &event.sftposgx);
    tree->Branch("sftposgy", &event.sftposgy);
    tree->Branch("sftposgz", &event.sftposgz);
    tree->Branch("sftpath", &event.sftpath);
    tree->Branch("sftbeta", &event.sftbeta);
    tree->Branch("sftmass", &event.sftmass);
  }

}
