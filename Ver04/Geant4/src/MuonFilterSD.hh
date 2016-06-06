/*
  MuonFilterSD.hh

  2016/2  K.Shirotori
*/

#ifndef MuonFilterSD_h
#define MuonFilterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MuonFilterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class MuonFilterSD : public G4VSensitiveDetector
{
public:
  MuonFilterSD( G4String name );
  ~MuonFilterSD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  MuonFilterHitsCollection *MuonFilterCollection;
};

#endif
