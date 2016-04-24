/*
  EventAction.cc

  2012/5  K.Shirotori
*/

#include "EventAction.hh"

EventAction::EventAction( Analysis *analysisManager, G4bool fAlwaysDraw )
  : G4UserEventAction(), anaMan(analysisManager), 
    fAlwaysDraw_(fAlwaysDraw)
{}

EventAction::~EventAction()
{}

