/*
  RICHSD.hh

  2015/12  K.Shirotori
*/

#ifndef RICHSD_h
#define RICHSD_h 1

#include "G4VSensitiveDetector.hh"
#include "RICHHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class RICHSD : public G4VSensitiveDetector
{
public:
  RICHSD( G4String name );
  ~RICHSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  RICHHitsCollection *RICHCollection;
};

#endif
