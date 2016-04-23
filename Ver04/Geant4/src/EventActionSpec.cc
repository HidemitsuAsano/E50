/*
  EventActionSpec.cc

  2012/5  K.Shirotori
*/

#include "EventActionSpec.hh"
#include "Analysis.hh"
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

#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4VTrajectory.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VTrajectoryPoint.hh"

#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>


EventActionSpec::
EventActionSpec( Analysis *analysisManager, G4bool fAlwaysDraw )
  : EventAction(analysisManager,fAlwaysDraw)
{}

EventActionSpec::~EventActionSpec()
{ 
}

void EventActionSpec::BeginOfEventAction( const G4Event *anEvent )
{
  G4int eventID = anEvent->GetEventID();
  if( eventID%1000==0 )
  //if( eventID%1000000==0 )
    G4cout << "Event: " << std::setw(8) << eventID << G4endl;
  
  G4SDManager *SDManager = G4SDManager::GetSDMpointer();
  if( TargetcolID<0    || 
      T0colID<0        ||
      BTrackercolID<0  || 
      SITrackercolID<0 ||
      AOTrackercolID<0 ||
      IITrackercolID<0 || 
      IOTrackercolID<0 ||
      SOTrackercolID<0 || 
      TofcolID<0       || 
      ITofcolID<0      || 
      PADcolID<0       || 
      RICHcolID<0      || 
      PID1colID<0      || 
      PID2colID<0      || 
      MuonFiltercolID<0|| 
      VDcolID<0 ){
    TargetcolID     = SDManager->GetCollectionID( G4String( "TargetCollection" ) );
    T0colID         = SDManager->GetCollectionID( G4String( "T0Collection" ) ); 
    BTrackercolID   = SDManager->GetCollectionID( G4String( "BTrackerCollection" ) );
    SITrackercolID  = SDManager->GetCollectionID( G4String( "SITrackerCollection" ) );
    AOTrackercolID  = SDManager->GetCollectionID( G4String( "AOTrackerCollection" ) );
    IITrackercolID  = SDManager->GetCollectionID( G4String( "IITrackerCollection" ) );
    IOTrackercolID  = SDManager->GetCollectionID( G4String( "IOTrackerCollection" ) );
    SOTrackercolID  = SDManager->GetCollectionID( G4String( "SOTrackerCollection" ) );
    TofcolID        = SDManager->GetCollectionID( G4String( "TofCollection" ) ); 
    ITofcolID       = SDManager->GetCollectionID( G4String( "ITofCollection" ) ); 
    PADcolID        = SDManager->GetCollectionID( G4String( "PADCollection" ) ); 
    RICHcolID       = SDManager->GetCollectionID( G4String( "RICHCollection" ) ); 
    PID1colID       = SDManager->GetCollectionID( G4String( "PID1Collection" ) ); 
    PID2colID       = SDManager->GetCollectionID( G4String( "PID2Collection" ) ); 
    MuonFiltercolID = SDManager->GetCollectionID( G4String( "MuonFilterCollection" ) ); 
    VDcolID         = SDManager->GetCollectionID( G4String( "VDCollection" ) ); 
  }

#if 0
  G4cout << "SpecEventActionSpec::BeginOfEventAction "
	 << "Collection IDs :" 
	 << TargercolID << " " 
	 << T0colID << " " 
	 << SITrackercolID << " " 
	 << AOTrackercolID << " " 
	 << IITrackercolID << " " 
	 << IOTrackercolID << " " 
	 << SOTrackercolID << " " 
	 << TofcolID << " " 
	 << ITofcolID << " " 
	 << PADcolID << " " 
	 << RICHcolID << " " 
	 << PID1colID << " " 
	 << PID2colID << " " 
	 << MuonFiltercolID << " " 
	 << VDcolID << G4endl;
#endif
  
  if( anaMan ) anaMan->BeginOfEvent( anEvent );
}

void EventActionSpec::EndOfEventAction( const G4Event *anEvent )
{
  G4VVisManager *visManager = G4VVisManager::GetConcreteInstance();

  G4bool fTriggered = false;
  if(anaMan){
    anaMan->EndOfEvent( anEvent );
    fTriggered = anaMan->GetTriggerStatus();
  }
  //if( visManager && ( fAlwaysDraw_ || fTriggered ) ){
  if( visManager ){
    //G4UImanager *UI = G4UImanager::GetUIpointer();
    // Trajectories
    G4TrajectoryContainer *trajectoryContainer =
      anEvent->GetTrajectoryContainer();
    if( trajectoryContainer ){
      G4int nTrajectories = trajectoryContainer->entries();
      for( G4int i=0; i<nTrajectories; i++ ){
	G4VTrajectory *trj = (*trajectoryContainer)[i];
	trj->DrawTrajectory(10);
#if 0
	G4int nPoint=trj->GetPointEntries();
	for( G4int j=0; j<nPoint; ++j ){
	  G4VTrajectoryPoint *tp=trj->GetPoint(j);
	  G4ThreeVector position=tp->GetPosition();
	  G4cout << position << G4endl;
	}
#endif
      }
    }
    // Hits 
//     G4HCofThisEvent *HCE = anEvent->GetHCofThisEvent();
//     if( HCE ){ 
//       if( TargetcolID>=0 ){
// 	G4VHitsCollection *TargetHC  = HCE->GetHC( TargetcolID );
// 	if( TargetHC )  TargetHC->DrawAllHits();
//       }
//       if( TrackercolID>=0 ){
// 	G4VHitsCollection *TrackerHC  = HCE->GetHC( TrackercolID );
// 	if( TrackerHC )  TrackerHC->DrawAllHits();
//       }
//       if( TofcolID>=0 ){
// 	G4VHitsCollection *TofHC = HCE->GetHC( TofcolID );
// 	if( TofHC ) TofHC->DrawAllHits();

//       }
//       if( T0colID>=0 ){
// 	G4VHitsCollection *T0HC  = HCE->GetHC( T0colID );
// 	if( T0HC )  T0HC->DrawAllHits();
//       }     
//       if( VDcolID>=0 ){
// 	G4VHitsCollection *VDHC  = HCE->GetHC( VDcolID );
// 	if( VDHC )  VDHC->DrawAllHits();
//       }     
    //     UI->ApplyCommand( "/vis/viewer/update" );
  }
}
