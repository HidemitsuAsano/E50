/*
  PID2SD.hh

  2016/2  K.Shirotori
*/

#ifndef PID2SD_h
#define PID2SD_h 1

#include "G4VSensitiveDetector.hh"
#include "PID2Hit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PID2SD : public G4VSensitiveDetector
{
public:
  PID2SD( G4String name );
  ~PID2SD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  PID2HitsCollection *PID2Collection;
};

#endif
