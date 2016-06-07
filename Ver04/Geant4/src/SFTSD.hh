/*
  SFTSD.hh

  2016/4  K.Shirotori
*/

#ifndef SFTSD_h
#define SFTSD_h 1

#include "G4VSensitiveDetector.hh"
#include "SFTHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class SFTSD : public G4VSensitiveDetector
{
public:
  SFTSD( G4String name );
  ~SFTSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();
  
private:
  SFTHitsCollection *SFTCollection;
};

#endif
