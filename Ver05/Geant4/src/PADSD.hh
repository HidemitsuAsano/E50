/*
  PADSD.hh

  2015/12  K.Shirotori
*/

#ifndef PADSD_h
#define PADSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PADHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class PADSD : public G4VSensitiveDetector
{
public:
  PADSD( G4String name );
  ~PADSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  PADHitsCollection *PADCollection;
};

#endif
