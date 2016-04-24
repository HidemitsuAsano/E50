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
  G4int T0colID;
  G4int SFTcolID;
};

#endif

  
