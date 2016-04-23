/*
  s_AFT.hh

  Simple AFT
  -> Just plate + Frame

  2015/12  K.Shirotori
*/

#ifndef s_AFT_h
#define s_AFT_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_AFT
{
public:  
  s_AFT( const G4String & Cname,
	 G4LogicalVolume *pMother, 
	 const G4RotationMatrix & rotMat, 
	 const G4ThreeVector & gPos1, 
	 const G4ThreeVector & gPos2,
	 const G4ThreeVector & OffsetLocal,
	 G4int id1,  
	 G4int id2, 
	 G4Material *matScin, 
	 G4Material *matFrame, 
	 G4Material *matArea, 
	 G4Material *matBox );
  ~s_AFT() {}; 

private:
  s_AFT( const s_AFT & );
  s_AFT & operator = ( const s_AFT & );
  
  G4String Cname_;
  G4LogicalVolume *logBox, *logFrame, *logArea;
  G4LogicalVolume *logLayerX;
  
public:
  void SetVisAttributes( const G4VisAttributes *attLayer,
			 const G4VisAttributes *attFrame,
			 const G4VisAttributes *attArea,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};
#endif
