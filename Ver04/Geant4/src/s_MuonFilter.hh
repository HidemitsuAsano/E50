/*
  s_MuonFilter.hh
  
  Simple MuonFilter -> Just plate

  2016/2  K.Shirotori
*/

#ifndef s_MuonFilter_h
#define s_MuonFilter_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class s_MuonFilter
{
public:
  s_MuonFilter( const G4String & Cname, 
		G4LogicalVolume *pMother,	    
		const G4RotationMatrix & rotMat, 
		const G4ThreeVector & gPosU, 
		const G4ThreeVector & gPosD, 
		const G4ThreeVector & OffSetLocal,
		G4Material *matScinti,
 		G4Material *matIron,
		G4Material *matBox );
  ~s_MuonFilter() {};
  
private:
  s_MuonFilter( const s_MuonFilter & );
  s_MuonFilter & operator = (const s_MuonFilter & );

  G4String Cname_;
  G4LogicalVolume *logBox, *logScinU, *logScinD;
  G4LogicalVolume *logIron;

public:
  void SetVisAttributes( const G4VisAttributes *attScin,
			 const G4VisAttributes *attIron,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );
};

#endif
