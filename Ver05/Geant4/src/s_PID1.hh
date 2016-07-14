/*
  s_PID1.hh
  
  Simple PID1 -> Just plate and box

  2015/12  K.Shirotori
*/

#ifndef s_PID1_h
#define s_PID1_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_PID1
{
public:
  s_PID1( const G4String & Cname, 
	  G4LogicalVolume *pMother,	    
	  const G4RotationMatrix & rotMat, 
	  const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	  G4Material *matAerogel, 
	  G4Material *matGas, 
	  G4Material *matFrame, 
	  G4Material *matBox );
  ~s_PID1() {};

private:
  s_PID1( const s_PID1 & );
  s_PID1 & operator = (const s_PID1 & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logSeg;

public:
  void SetVisAttributes( const G4VisAttributes *attSeg,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
