/*
  r_SFT_D.hh

  Segmented SFT
  -> Fiber + Frame

  2016/4  K.Shirotori
*/

#ifndef r_SFT_D_h
#define r_SFT_D_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class r_SFT_D
{
public:  
  r_SFT_D( const G4String & Cname,
	   G4VPhysicalVolume *pMother, 
	   const G4RotationMatrix & rotMat, 
	   const G4ThreeVector *gPos1, 
	   const G4ThreeVector & OffsetLocal,
	   G4int *id, 
     const G4double *LzLayer,
	   G4Material *matScin, 
	   G4Material *matPMMA, 
	   G4Material *matFP, 
	   G4Material *matFrame, 
	   G4Material *matArea, 
	   G4Material *matBox );
  ~r_SFT_D() {}; 
     
private:
  r_SFT_D( const r_SFT_D & );
  r_SFT_D & operator = ( const r_SFT_D & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logLayer_core;
  G4LogicalVolume *logLayer_innerclad;
  G4LogicalVolume *logLayer_outerclad;
  
public:
  void SetVisAttributes( const G4VisAttributes *attLayer,
      const G4VisAttributes *attInnerClad,
      const G4VisAttributes *attOuterClad,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );



};
#endif
