/*
  s_T0Wall.hh
  
  Simple T0 -> Just plate

  2015/12  K.Shirotori
*/

#ifndef s_T0Wall_h
#define s_T0Wall_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_T0Wall
{
public:
  s_T0Wall( const G4String & Cname, 
	    G4VPhysicalVolume *pMother,
	    const G4RotationMatrix & rotMat, 
	    const G4ThreeVector & gPos, 
	    const G4ThreeVector & OffSetLocal,
	    G4Material *matScinti,
	    G4Material *matBox );
  ~s_T0Wall() {};

private:
  s_T0Wall( const s_T0Wall & );
  s_T0Wall & operator = (const s_T0Wall & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logScin;

public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
