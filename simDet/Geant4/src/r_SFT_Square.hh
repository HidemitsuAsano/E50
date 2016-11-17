/*
  r_SFT_Square.hh

  Segmented SFT
  -> Fiber + Frame

  2016/4  K.Shirotori
*/

#ifndef r_SFT_Square_h
#define r_SFT_Square_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class r_SFT_Square
{
public:  
  r_SFT_Square( const G4String & Cname,
	   G4VPhysicalVolume *pMother, 
	   const G4RotationMatrix & rotMat, 
	   const G4ThreeVector *gPos1, 
	   const G4ThreeVector & OffsetLocal,
	   const G4int *detid, 
	   const G4Material *matScin, 
	   const G4Material *matPMMA, 
	   const G4Material *matFrame, 
	   const G4Material *matArea, 
	   const G4Material *matBox );
  ~r_SFT_Square() {}; 
  
private:
  r_SFT_Square( const r_SFT_Square & );
  r_SFT_Square & operator = ( const r_SFT_Square & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logLayerX_core, *logLayerU_core, *logLayerV_core;
  G4LogicalVolume *logLayerX_clad, *logLayerU_clad, *logLayerV_clad;
  
public:
  void SetVisAttributes( const G4VisAttributes *attLayer,
      const G4VisAttributes *attClad,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};
#endif
