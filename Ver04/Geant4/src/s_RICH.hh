/*
  s_RICH.hh
  
  Simple RICH -> Just plate and box

  2015/12  K.Shirotori
*/

#ifndef s_RICH_h
#define s_RICH_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_RICH
{
public:
  s_RICH( const G4String & Cname, 
	  G4LogicalVolume *pMother,	    
	  const G4RotationMatrix & rotMat, 
	  const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	  G4Material *matAerogel, 
	  G4Material *matGas, 
	  G4Material *matFrame, 
	  G4Material *matBox );
  ~s_RICH() {};

private:
  s_RICH( const s_RICH & );
  s_RICH & operator = (const s_RICH & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logGas, *logAerogel;

public:
  void SetVisAttributes( const G4VisAttributes *attAerogel,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
