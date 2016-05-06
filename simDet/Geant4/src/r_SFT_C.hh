/*
  r_SFT_C.hh

  Segmented SFT 0.5mm
  -> Fiber + Frame

  May 2016  H.Asano
*/

#ifndef r_SFT_C_h
#define r_SFT_C_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class r_SFT_C
{
public:  
  r_SFT_C( const G4String & Cname,
	   G4VPhysicalVolume *pMother, 
	   const G4RotationMatrix & rotMat, 
	   const G4ThreeVector *gPos, 
	   const G4ThreeVector & OffsetLocal,
	   G4int *detid, 
	   G4Material *matScin, 
	   G4Material *matPMMA, 
	   G4Material *matFrame, 
	   G4Material *matArea, 
	   G4Material *matBox );
  ~r_SFT_C() {}; 
  
private:
  r_SFT_C( const r_SFT_C & );
  r_SFT_C & operator = ( const r_SFT_C & );
  
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
