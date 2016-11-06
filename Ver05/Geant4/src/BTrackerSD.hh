/*
  BTrackerSD.hh

  2013/1  K.Shirotori
*/

#ifndef BTrackerSD_h
#define BTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "BTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class BTrackerSD : public G4VSensitiveDetector
{
public:
  BTrackerSD( G4String name );
  ~BTrackerSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  BTrackerHitsCollection *BTrackerCollection;
};

#endif
