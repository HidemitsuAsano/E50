/*
  SteppingAction.cc

  2016/1  K.Shirotori
*/

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"

#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "ConfMan.hh"

SteppingAction::SteppingAction( DetectorConstruction* det,
				EventAction* evt )
  : detector(det), eventaction(evt)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* theStep )
{ 
  ConfMan *confMan = ConfMan::GetConfManager();
  int StepFlag = confMan->StepFlag();

  if(!StepFlag) return;

  G4Track * theTrack = theStep->GetTrack();

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(theStep->GetPreStepPoint()->GetTouchable());
    
  G4VPhysicalVolume* physVol = theTouchable->GetVolume(); 

  G4ThreeVector pos = theTrack->GetPosition();
  G4StepPoint *preStepPoint = theStep->GetPreStepPoint();
  G4TouchableHandle Touchable = preStepPoint->GetTouchableHandle();
  G4ThreeVector posl = Touchable->GetHistory()->
    GetTopTransform().TransformPoint( pos );

  // All particle which hits Spec is killed.
  if( GetFStopSPEC(StepFlag) && 
      detector->IsVolumeStopper(physVol) ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  G4String particle = theTrack->GetDefinition()->GetParticleName();

  // Gamma ray which hits Spec is killed.
  if( GetFStopSPECGam(StepFlag) && particle=="gamma" &&
      detector->IsVolumeStopper(physVol) ){
    theTrack->SetTrackStatus( fStopAndKill );
    return;
  }

  // All neutrino are killed.
  if( GetFStopNu(StepFlag) &&
      ( particle =="anti_nu_e" || particle =="anti_nu_mu" 
	|| particle =="nu_e" || particle =="nu_mu" ) ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  // Gamma is killed.
  if( GetFStopGam(StepFlag) && particle =="gamma" ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }
  // e+/e- are killed.
  if( GetFStopE(StepFlag) && ( particle =="e-" || particle =="e+" ) ){
    theTrack->SetTrackStatus(fStopAndKill);
    return;
  }

  if( theTrack->GetTrackLength() >= 10.0*m ){
    theTrack->SetTrackStatus(fStopAndKill);      
    return;
  }

  // // Kaon+ is killed.
  // if( GetFStopGam(StepFlag) && particle =="kaonP1" ){
  //   theTrack->SetTrackStatus(fStopAndKill);
  //   return;
  // }

  // // Pion- is killed.
  // if( GetFStopGam(StepFlag) && particle =="pionM1" ){
  //   theTrack->SetTrackStatus(fStopAndKill);
  //   return;
  // }

//   if( theTrack->GetStepLength() < 0.0000001*mm ){
//     theTrack->SetTrackStatus(fStopAndKill);      
//     return;
//   }

  return;
}

