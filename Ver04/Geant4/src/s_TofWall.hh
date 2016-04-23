/*
  s_TofWall.hh
  
  Simple TOF -> Just plate

  2015/12  K.Shirotori
*/

#ifndef s_TofWall_h
#define s_TofWall_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_TofWall
{
public:
  s_TofWall( const G4String & Cname, 
	     G4LogicalVolume *pMother,	    
	     const G4RotationMatrix & rotMat, 
	     const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	     G4Material *matScinti, 
	     G4Material *matBox );
  ~s_TofWall() {};

private:
  s_TofWall( const s_TofWall & );
  s_TofWall & operator = (const s_TofWall & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logScin;

public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
