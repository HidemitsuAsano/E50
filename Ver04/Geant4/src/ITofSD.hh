/*
  ITofSD.hh

  2012/8  K.Shirotori
*/

#ifndef ITofSD_h
#define ITofSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ITofHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ITofSD : public G4VSensitiveDetector
{
public:
  ITofSD( G4String name );
  ~ITofSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  ITofHitsCollection *ITofCollection;
};

#endif
