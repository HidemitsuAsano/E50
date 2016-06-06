/*
  VDSD.hh

  2012/5  K.Shirotori
*/

#ifndef VDSD_h
#define VDSD_h 1

#include "G4VSensitiveDetector.hh"
#include "VDHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class VDSD : public G4VSensitiveDetector
{
public:
  VDSD( G4String name );
  ~VDSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  VDHitsCollection *VDCollection;
};

#endif
