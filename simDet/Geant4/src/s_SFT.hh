/*
  s_SFT.hh

  Simple SFT
  -> Just plate + Frame

  2015/12  K.Shirotori
  
  
  2016/5  updated by H.Asano
*/

#ifndef s_SFT_h
#define s_SFT_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_SFT
{
public:  
  s_SFT( const G4String & Cname,
	 G4VPhysicalVolume *pMother, 
	 const G4RotationMatrix & rotMat, 
	 const G4ThreeVector *gPos, 
	 const G4ThreeVector & OffsetLocal,
	 const G4int *detid, 
	 G4Material *matScin, 
	 G4Material *matFrame, 
	 G4Material *matArea, 
	 G4Material *matBox );
  ~s_SFT() {}; 

private:
  s_SFT( const s_SFT & );
  s_SFT & operator = ( const s_SFT & );
  
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
