/*
  AOTrackerSD.hh

  2016/2  K.Shirotori
*/

#ifndef AOTrackerSD_h
#define AOTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "AOTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class AOTrackerSD : public G4VSensitiveDetector
{
public:
  AOTrackerSD( G4String name );
  ~AOTrackerSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  AOTrackerHitsCollection *AOTrackerCollection;
};

#endif
