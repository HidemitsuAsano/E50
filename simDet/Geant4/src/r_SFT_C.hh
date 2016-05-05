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
	   const G4ThreeVector & gPos1, 
	   const G4ThreeVector & gPos2,
	   const G4ThreeVector & gPos3, 
	   const G4ThreeVector & gPos4, 
	   const G4ThreeVector & gPos5, 
	   const G4ThreeVector & gPos6, 
	   const G4ThreeVector & gPos7, 
	   const G4ThreeVector & gPos8, 
	   const G4ThreeVector & gPos9, 
	   const G4ThreeVector & gPos10, 
	   const G4ThreeVector & gPos11, 
	   const G4ThreeVector & gPos12, 
	   const G4ThreeVector & OffsetLocal,
	   G4int id1,  G4int id2,  G4int id3, 
	   G4int id4,  G4int id5,  G4int id6,
	   G4int id7,  G4int id8,  G4int id9,
	   G4int id10, G4int id11, G4int id12,
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
