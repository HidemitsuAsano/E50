/*
  SOTrackerSD.hh

  2013/1  K.Shirotori
*/

#ifndef SOTrackerSD_h
#define SOTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SOTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SOTrackerSD : public G4VSensitiveDetector
{
public:
  SOTrackerSD( G4String name );
  ~SOTrackerSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  SOTrackerHitsCollection *SOTrackerCollection;
};

#endif
