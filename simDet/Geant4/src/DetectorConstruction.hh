/*
  DetectorConstruction.hh

  2016/4  K.Shirotori
*/

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"

#include <string>

class G4Material;
class MaterialList;

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

protected:
  MaterialList *mList_;

  MaterialList *DefineMaterials( void );
  virtual G4VPhysicalVolume *ConstructPayload( void ) = 0;
};

#endif
