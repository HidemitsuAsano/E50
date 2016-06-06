/*
  SOTrackerSD.cc

  2013/1  K.Shirotori
*/

#include "SOTrackerSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4TouchableHandle.hh"

#include "ConfMan.hh"

const G4double DEThreshold = 0.001*keV;
const double PositionSeparationThreshold = 0.1*mm;
const double TimeSeparationThreshold     = 5.0*ns;

SOTrackerSD::SOTrackerSD( G4String name )
  : G4VSensitiveDetector(name)
{
  collectionName.insert( G4String( "SOTrackerCollection" ) );
}

SOTrackerSD::~SOTrackerSD()
{
}

void SOTrackerSD::Initialize( G4HCofThisEvent *HCE )
{
  static int HCID = -1;
  SOTrackerCollection =
    new SOTrackerHitsCollection( SensitiveDetectorName,
			       collectionName[0] );
  if( HCID<0 )
    HCID = GetCollectionID(0);
  
  HCE->AddHitsCollection( HCID, SOTrackerCollection );
}

G4bool SOTrackerSD::ProcessHits( G4Step *aStep,
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
  G4int hitLayer = theTouchable->GetReplicaNumber();
  G4int hitSegment = vol->GetCopyNo();

  G4int nHits = SOTrackerCollection->entries();
  G4ThreeVector hitpos = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime = aTrack->GetGlobalTime();
  G4int trackNo = aTrack->GetTrackID();
  G4ThreeVector hitmom = aTrack->GetMomentum();
  G4double path = aTrack->GetTrackLength();
  G4double beta = aStep->GetPreStepPoint()->GetBeta(); 
  G4double hitmass = aTrack->GetDefinition()->GetPDGMass();
  
  G4String PartName= aTrack->GetDefinition()->GetParticleName();
  G4ThreeVector hitposl = theTouchable->GetHistory()->
    GetTopTransform().TransformPoint( hitpos );
  G4ThreeVector hitmoml = theTouchable->GetHistory()->
    GetTopTransform().TransformAxis( hitmom );
  
  for( G4int i=0; i<nHits; ++i ){
    SOTrackerHit *aHit = (*SOTrackerCollection)[i];
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
  G4int layer=0;
  for( G4int i=0; i<nHits; ++i ){
    SOTrackerHit *aHit=(*SOTrackerCollection)[i];
    id = aHit->GetTrackNo();
    layer = aHit->GetLayerID();
    if( id == trackNo && layer == hitLayer ) return true;
  }

  SOTrackerHit *aHit = new SOTrackerHit();
  aHit->SetPass();
  aHit->SetLayerID( hitLayer );
  aHit->SetSegmentID( hitSegment );
  aHit->SetEdep( edep );
  aHit->SetTime( hittime );
  aHit->SetPos( hitpos );
  aHit->SetMom( hitmom );
  aHit->SetLMom( hitmoml );
  aHit->SetTrackNo( trackNo );
  aHit->SetLocalPos( hitposl.x(), hitposl.z() );
  aHit->SetHitParticleName( PartName );
  aHit->SetPathLength( path );
  aHit->SetBeta( beta );
  aHit->SetHitParticleMass( hitmass );
  SOTrackerCollection->insert( aHit );

#if 0
  G4cout << "[SOTrackerSD] " 
	 << "Layer="    << hitLayer 
	 << " TrackID=" << trackNo 
	 << " Name=" << PartName 
	 << " Mass=" << hitmass //<< G4endl;
// 	 << " edep="    << edep/keV << "keV"  
// 	 << " G: "      << hitpos 
	 << " L: "      << hitposl 
// 	 << " P: "      << hitmom 
	 << G4endl;
#endif

  return true;
}

void SOTrackerSD::EndOfEvent( G4HCofThisEvent *HCE )
{
}

void SOTrackerSD::clear()
{
  G4int nHits = SOTrackerCollection->entries();
  for( G4int i=nHits-1; i>=0; --i )
    delete (*SOTrackerCollection)[i];
}

void SOTrackerSD::DrawAll() const
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();

  if( pVisManager ){
    G4int nHits = SOTrackerCollection->entries();
    for( G4int i=0; i<nHits; ++i )
      (*SOTrackerCollection)[i]->Draw();
  }
}

void SOTrackerSD::PrintAll() const
{
  G4int nHits = SOTrackerCollection->entries();
  for( G4int i=0; i<nHits; ++i)
    (*SOTrackerCollection)[i]->Print();
}

