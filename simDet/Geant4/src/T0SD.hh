/*
  T0SD.hh

  2012/5  K.Shirotori
*/

#ifndef T0SD_h
#define T0SD_h 1

#include "G4VSensitiveDetector.hh"
#include "T0Hit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class T0SD : public G4VSensitiveDetector
{
public:
  T0SD( G4String name );
  ~T0SD();

  void Initialize( G4HCofThisEvent *HCE );
  G4bool ProcessHits( G4Step *aStep, G4TouchableHistory *ROhist );
  void EndOfEvent( G4HCofThisEvent *HCE );

  void DrawAll() const;
  void PrintAll() const;
  void clear();

private:
  T0HitsCollection *T0Collection;
  G4int Verbosity_;
};

#endif
