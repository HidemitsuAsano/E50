/*
  IITrackerSD.hh

  2013/4  K.Shirotori
*/

#ifndef IITrackerSD_h
#define IITrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "IITrackerHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class IITrackerSD : public G4VSensitiveDetector
{
public:
  IITrackerSD( G4String name );
  ~IITrackerSD();
  
  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );
  
  void DrawAll() const;
  void PrintAll() const;
  void clear();
  
private:
  IITrackerHitsCollection *IITrackerCollection;
};

#endif
