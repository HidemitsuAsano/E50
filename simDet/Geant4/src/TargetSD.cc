/*
  TargetSD.cc

  2012/5  K.Shirotori
*/

#include "TargetSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"

#include "ConfMan.hh"

const double PositionSeparationThreshold = 2.0*cm;
const double TimeSeparationThreshold     = 5.0*ns;

TargetSD::TargetSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( G4String( "TargetCollection" ) );
}

TargetSD::~TargetSD()
{
}

void TargetSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  TargetCollection =
    new TargetHitsCollection( SensitiveDetectorName,
			     collectionName[0] );
  if( HCID<0 )
    HCID = GetCollectionID(0);
  
  HCE->AddHitsCollection( HCID, TargetCollection );
}

G4bool TargetSD::ProcessHits( G4Step *aStep,
			      G4TouchableHistory *ROhist )
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  ConfMan *confMan = ConfMan::GetConfManager();
  bool detEvt = confMan->DetectEveryThing();
  if( (!detEvt) && edep==0.0 ) return true;

  G4Track *aTrack = aStep->GetTrack();
  const G4VTouchable *theTouchable =
    aStep->GetPreStepPoint()->GetTouchable();
  G4VPhysicalVolume *vol=theTouchable->GetVolume();
  G4String hitName = vol->GetName();
  G4int hitLayer=0;
  if( hitName=="Target" )            hitLayer=1;
  else if( hitName=="TargetBox" )    hitLayer=2;
  G4int hitSegment = vol->GetCopyNo();

  G4int nHits = TargetCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime = aTrack->GetGlobalTime();
  G4int trackNo = aTrack->GetTrackID();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4double path = aTrack->GetTrackLength();
  G4double beta = aStep->GetPreStepPoint()->GetBeta(); 
  G4double hitmass = aTrack->GetDefinition()->GetPDGMass();

  G4ThreeVector hitmom2 = aStep->GetPreStepPoint()->GetMomentum();
  
  G4String PartName= aTrack->GetDefinition()->GetParticleName();
  G4ThreeVector hitposl = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint( hitpos );
  G4ThreeVector hitmoml = theTouchable->GetHistory()->
    GetTopTransform().TransformAxis( hitmom );

  for( G4int i=0; i<nHits; ++i ){
    TargetHit *aHit = (*TargetCollection)[i];
    if( hitLayer==aHit->GetLayerID() && 
	hitSegment==aHit->GetSegmentID() ){
      G4double time = aHit->GetTime();
      G4double id = aHit->GetTrackNo();
      if( fabs(hittime-time)<=TimeSeparationThreshold &&
	  id == trackNo )
	aHit->AddEdep( edep );
    }
  }

  G4int id=0;
  for( G4int i=0; i<nHits; ++i ){
    TargetHit *aHit=(*TargetCollection)[i];
    id = aHit->GetTrackNo();
    if( id == trackNo ) return true;
  }

  TargetHit *aHit = new TargetHit();
  aHit->SetPass();
  aHit->SetLayerID( hitLayer );
  aHit->SetSegmentID( hitSegment );
  aHit->SetEdep( edep );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom2 );
  //aHit->SetMom( hitmom );
  aHit->SetLMom( hitmoml );
  aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.x(), hitposl.y() );
  aHit->SetHitParticleName( PartName );
  aHit->SetPathLength( path );
  aHit->SetBeta( beta );
  aHit->SetHitParticleMass( hitmass );
  TargetCollection->insert( aHit );

 

#if 0
  G4cout << "[TargetSD] " 
	 << "Layer="    << hitLayer 
	 << " TrackID=" << trackNo 
	 << " Name=" << PartName 
	 << " Mass=" << hitmass << G4endl;
// 	 << " edep="    << edep/keV << "keV"  
// 	 << " G: "      << hitpos 
// 	 << " L: "      << hitposl 
// 	 << " P: "      << hitmom << G4endl;
#endif

  return true;
}

void TargetSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void TargetSD::clear()
{
  G4int nHits = TargetCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    delete (*TargetCollection)[i];
}

void TargetSD::DrawAll() const
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

  if( pVisManager ){
    G4int nHits = TargetCollection->entries();
    for( G4int i=0; i<nHits; ++i )
      (*TargetCollection)[i]->Draw();
  }
}

void TargetSD::PrintAll() const
{
  G4int nHits = TargetCollection->entries();
  for( G4int i=0; i<nHits; ++i)
    (*TargetCollection)[i]->Print();
}

