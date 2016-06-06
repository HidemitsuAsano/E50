/*
  r_T0Wall.hh

  Real T0

  2015/12  K.Shirotori
*/

#ifndef r_T0Wall_h
#define r_T0Wall_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class r_T0Wall
{
public:
  r_T0Wall( const G4String & Cname, 
	    G4LogicalVolume *pMother,
	    const G4RotationMatrix & rotMat, 
	    const G4ThreeVector & gPos, 
	    const G4ThreeVector & OffSetLocal,
	    G4Material *matScinti, 
	    G4Material *matMPPC,
	    G4Material *matPCB, 
	    G4Material *matBox );
  ~r_T0Wall() {};
private:
  r_T0Wall( const r_T0Wall & );
  r_T0Wall & operator = (const r_T0Wall & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logScin, *logMPPC, *logPCB;
public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attMPPC,
			 const G4VisAttributes *attPCM,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
