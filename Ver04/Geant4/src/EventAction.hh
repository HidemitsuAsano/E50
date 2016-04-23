/*
  EventAction.hh

  2012/5  K.Shirotori
*/

#ifndef EventAction_h 
#define EventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"

class Analysis;

class EventAction : public G4UserEventAction
{
public:
  EventAction( Analysis *analysisManager=0, G4bool fAlwaysDraw=false );
  ~EventAction();

private:
  EventAction( const EventAction & );
  EventAction & operator = ( const EventAction & );

public:
  virtual void BeginOfEventAction( const G4Event *anEvent )=0;
  virtual void EndOfEventAction( const G4Event *anEvent )=0;

  G4bool SetDrawFlag( G4bool NewVal ) 
  { G4bool OldVal=fAlwaysDraw_; fAlwaysDraw_=NewVal; return OldVal; }

protected:
  Analysis *anaMan;
  G4bool fAlwaysDraw_;

};

#endif
