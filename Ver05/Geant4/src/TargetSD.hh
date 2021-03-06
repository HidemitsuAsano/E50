/*
  TargetSD.hh

  2012/5  K.Shirotori
*/

#ifndef TargetSD_h
#define TargetSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TargetHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class TargetSD : public G4VSensitiveDetector
{
public:
  TargetSD( G4String name );
  ~TargetSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  TargetHitsCollection *TargetCollection;
};

#endif
