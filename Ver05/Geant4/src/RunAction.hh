/*
  RunAction.hh

  2012/5  K.Shirotori
*/

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

class G4Run;
class Analysis;

class RunAction : public G4UserRunAction
{
public:
  RunAction( Analysis *analysisManager=0 );
  ~RunAction();

public:
  void BeginOfRunAction( const G4Run *aRun );
  void EndOfRunAction( const G4Run *aRun );

private:
  Analysis *anaMan;
};

#endif
