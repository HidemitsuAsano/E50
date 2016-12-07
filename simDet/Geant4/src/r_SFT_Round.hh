/*
  r_SFT_Round.hh

  Segmented SFT
  -> Fiber + Frame

  2016/4  K.Shirotori
  
  2016/5 updated by H. Asano
*/

#ifndef r_SFT_Round_h
#define r_SFT_Round_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class r_SFT_Round
{
public:  
  r_SFT_Round( const G4String & Cname,
	   G4VPhysicalVolume *pMother, 
	   const G4RotationMatrix & rotMat, 
	   const G4ThreeVector *gPos, 
	   const G4ThreeVector & OffsetLocal,
	   const G4int *detid, 
     const G4double *LzLayer,
	   G4Material *matScin, 
	   G4Material *matPMMA, 
	   G4Material *matFP, 
	   G4Material *matFrame, 
	   G4Material *matArea, 
	   G4Material *matBox );
  ~r_SFT_Round() {}; 
     
private:
  r_SFT_Round( const r_SFT_Round & );
  r_SFT_Round & operator = ( const r_SFT_Round & );
  G4String Cname_;
  G4int Verbosity_; 
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
