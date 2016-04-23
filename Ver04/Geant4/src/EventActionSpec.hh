/*
  EventActionSpec.hh

  2012/5  K.Shirotori
*/

#ifndef EventActionSpec_h 
#define EventActionSpec_h 1

#include "EventAction.hh"

class EventActionSpec : public EventAction
{
public:
  EventActionSpec( Analysis *analysisManager=0, 
		   G4bool fAlwaysDraw=false );
  ~EventActionSpec();
  
public:
  void BeginOfEventAction( const G4Event *anEvent );
  void EndOfEventAction( const G4Event *anEvent );
  
private:
  G4int TargetcolID;
  G4int T0colID;
  G4int BTrackercolID;
  G4int SITrackercolID;
  G4int AOTrackercolID;
  G4int IITrackercolID;
  G4int IOTrackercolID;
  G4int SOTrackercolID;
  G4int TofcolID;
  G4int ITofcolID;
  G4int PADcolID;
  G4int RICHcolID;
  G4int PID1colID;
  G4int PID2colID;
  G4int MuonFiltercolID;
  G4int VDcolID;
};

#endif

  
