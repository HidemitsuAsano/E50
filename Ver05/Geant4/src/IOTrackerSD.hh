/*
  IOTrackerSD.hh

  2013/4  K.Shirotori
*/

#ifndef IOTrackerSD_h
#define IOTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "IOTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class IOTrackerSD : public G4VSensitiveDetector
{
public:
  IOTrackerSD( G4String name );
  ~IOTrackerSD();
  
  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );
  
  void DrawAll() const;
  void PrintAll() const;
  void clear();
  
private:
  IOTrackerHitsCollection *IOTrackerCollection;
};

#endif
