/*
  s_ITofWall.hh
  
  Simple ITOF -> Just plate

  2015/12  K.Shirotori
*/

#ifndef s_ITofWall_h
#define s_ITofWall_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_ITofWall
{
public:
  s_ITofWall( const G4String & Cname, 
	      G4LogicalVolume *pMother,	    
	      const G4RotationMatrix & rotMat, 
	      const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	      G4Material *matScinti, 
	      G4Material *matBox );
  ~s_ITofWall() {};

private:
  s_ITofWall( const s_ITofWall & );
  s_ITofWall & operator = (const s_ITofWall & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logScin;

public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
