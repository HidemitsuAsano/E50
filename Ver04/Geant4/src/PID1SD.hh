/*
  PID1SD.hh

  2015/12  K.Shirotori
*/

#ifndef PID1SD_h
#define PID1SD_h 1

#include "G4VSensitiveDetector.hh"
#include "PID1Hit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PID1SD : public G4VSensitiveDetector
{
public:
  PID1SD( G4String name );
  ~PID1SD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  PID1HitsCollection *PID1Collection;
};

#endif
