/*
  LiqTarget.hh

  2012/5  K.Shirotori
*/

#ifndef LiqTarget_h
#define LiqTarget_h 1

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4Material;
class G4VisAttributes;
class G4VSensitiveDetector;
class G4LogicalVolume;
class G4VPhysicalVolume;

class LiqTarget
{
public:
  LiqTarget( const G4String & Cname, G4VPhysicalVolume *pMother,
	     const G4RotationMatrix & rotMat, 
	     const G4ThreeVector & gPos, const G4ThreeVector & OffSetLocal,
	     G4Material *matTarget, G4Material *matTargCell,
	     G4Material *matVessel, G4Material *matWindow,
	     G4Material *matVacuum, G4Material *matBox );
  ~LiqTarget() {};
private:
  LiqTarget( const LiqTarget & );
  LiqTarget & operator = ( const LiqTarget & );

  G4String Cname_;
  G4LogicalVolume *logTarg, *logCell, *logVessel, *logVacuum,
    *logWindowIn, *logWindowOut, *logBox;
  double length_, diameter_;
public:
  void SetVisAttributes( const G4VisAttributes *attTarg,
			 const G4VisAttributes *attCell,
			 const G4VisAttributes *attVessel,
			 const G4VisAttributes *attVacuum,
			 const G4VisAttributes *attWindow,
			 const G4VisAttributes *attBox );
  void SetSensitiveDetector( G4VSensitiveDetector *SD );

  double Length() const { return length_; }
  double PhiSize() const { return diameter_; }
};

#endif 
