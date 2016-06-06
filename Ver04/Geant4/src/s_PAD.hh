/*
  s_PAD.hh
  
  Simple PAD -> Just plate

  2015/12  K.Shirotori
*/

#ifndef s_PAD_h
#define s_PAD_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_PAD
{
public:
  s_PAD( const G4String & Cname, 
	 G4LogicalVolume *pMother,	    
	 const G4RotationMatrix & rotMat, 
	 const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	 G4Material *matScinti, 
	 G4Material *matBox );
  ~s_PAD() {};

private:
  s_PAD( const s_PAD & );
  s_PAD & operator = (const s_PAD & );

  G4String Cname_; 
  G4LogicalVolume *logBox, *logScin;

public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
