/*
  s_ITWC1.hh

  Simple Iinternal WC 1
  -> Just plate + Frame

  2015/12  K.Shirotori
*/

#ifndef s_ITWC1_h
#define s_ITWC1_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_ITWC1
{
public:  
  s_ITWC1( const G4String & Cname,
	   G4LogicalVolume *pMother, 
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
	   G4Material *matWire, 
	   G4Material *matFrame, 
	   G4Material *matGas, 
	   G4Material *matBox );
  ~s_ITWC1() {}; 

private:
  s_ITWC1( const s_ITWC1 & );
  s_ITWC1 & operator = ( const s_ITWC1 & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logArea, *logFrame;
  G4LogicalVolume *logLayerX, *logLayerU, *logLayerV;
  
public:
  void SetVisAttributes( const G4VisAttributes *attLayer,
			 const G4VisAttributes *attFrame,
                         const G4VisAttributes *attArea,
                         const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};
#endif
