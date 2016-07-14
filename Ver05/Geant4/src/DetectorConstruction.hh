/*
  DetectorConstruction.hh

  2012/5  K.Shirotori
*/

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include <string>

class G4Material;
class MaterialList;
class G4MagneticField;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

private:
  DetectorConstruction( const DetectorConstruction & );
  DetectorConstruction & operator = ( const DetectorConstruction & );

public:
  virtual G4VPhysicalVolume *Construct();
  virtual G4ThreeVector TargetPosition( void ) const = 0;
  G4double TargetLength( void ) const { return TargetSizeZ(); }
  virtual G4double TargetAngle( void ) const = 0;
  virtual G4double TargetSizeX( void ) const = 0;
  virtual G4double TargetSizeY( void ) const = 0;
  virtual G4double TargetSizeZ( void ) const = 0;
  virtual G4bool IsVolumeStopper( G4VPhysicalVolume *physVol ) const;

protected:
  MaterialList *mList_;
  G4MagneticField *field_;

  MaterialList *DefineMaterials( void );
  virtual G4VPhysicalVolume *ConstructPayload( void ) = 0;
  virtual void MakeSpecMagnet( G4VPhysicalVolume *pMother,
			       G4Material *matGap );
  G4MagneticField * MakeUniformMagField( G4double Bz );
  G4MagneticField * MakeMagFieldFromMap( const std::string &filename,
					 G4double NormFac=1.0 );
};

#endif
