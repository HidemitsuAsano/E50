/*
  EventActionSpec.cc

  2012/5  K.Shirotori
*/

#include "EventActionSpec.hh"
#include "Analysis.hh"
#include "T0Hit.hh"
#include "SFTHit.hh"

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
  if( T0colID<0        ||
      SFTcolID<0 ){
    T0colID   = SDManager->GetCollectionID( G4String( "T0Collection" ) ); 
    SFTcolID  = SDManager->GetCollectionID( G4String( "SFTCollection" ) );
  }

#if 0
  G4cout << "SpecEventActionSpec::BeginOfEventAction "
	 << "Collection IDs :" 
	 << T0colID << " " 
	 << SFTcolID << " " 
	 << G4endl;
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
