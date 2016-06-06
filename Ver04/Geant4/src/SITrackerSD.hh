/*
  SITrackerSD.hh

  2013/1  K.Shirotori
*/

#ifndef SITrackerSD_h
#define SITrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SITrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SITrackerSD : public G4VSensitiveDetector
{
public:
  SITrackerSD( G4String name );
  ~SITrackerSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  SITrackerHitsCollection *SITrackerCollection;
};

#endif
