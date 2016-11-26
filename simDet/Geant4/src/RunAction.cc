/*
  RunAction.cc

  2012/5  K.Shirotori
*/

#include "RunAction.hh"
#include "Analysis.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"

RunAction::RunAction( Analysis *analysisManager )
  : G4UserRunAction(), anaMan(analysisManager)
{}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction( const G4Run *aRun )
{
  if( G4VVisManager::GetConcreteInstance() ){
   // G4UImanager *UI = G4UImanager::GetUIpointer();
    //    UI->ApplyCommand( "/vis/scene/notifyHandlers" );
  }
  if( anaMan ) anaMan->BeginOfRun( aRun );
}

void RunAction::EndOfRunAction( const G4Run *aRun )
{
  if( G4VVisManager::GetConcreteInstance() ){
  //  G4UImanager *UI = G4UImanager::GetUIpointer();
    //    UI->ApplyCommand( "/vis/viewer/update" );
  }

  if( anaMan ) anaMan->EndOfRun( aRun );
} 

