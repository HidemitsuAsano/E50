/*
  s_BFT.hh

  Simple BFT
  -> Just plate + Frame

  2015/12  K.Shirotori
*/

#ifndef s_BFT_h
#define s_BFT_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_BFT
{
public:  
  s_BFT( const G4String & Cname,
	 G4LogicalVolume *pMother, 
	 const G4RotationMatrix & rotMat, 
	 const G4ThreeVector & gPos1, 
	 const G4ThreeVector & gPos2,
	 const G4ThreeVector & gPos3, 
	 const G4ThreeVector & gPos4, 
	 const G4ThreeVector & gPos5, 
	 const G4ThreeVector & gPos6, 
	 const G4ThreeVector & OffsetLocal,
	 G4int id1,  G4int id2,  G4int id3, 
	 G4int id4,  G4int id5,  G4int id6,
	 G4Material *matScin, 
	 G4Material *matFrame, 
	 G4Material *matArea, 
	 G4Material *matBox );
  ~s_BFT() {}; 

private:
  s_BFT( const s_BFT & );
  s_BFT & operator = ( const s_BFT & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logLayerX, *logLayerU, *logLayerV;
  
public:
  void SetVisAttributes( const G4VisAttributes *attLayer,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};
#endif
