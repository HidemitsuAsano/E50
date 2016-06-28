/*
  DetectorConstructionSpec.hh

  2012/5  K.Shirotori
*/

#ifndef DetectorConstructionSpec_h
#define DetectorConstructionSpec_h 1

#include "DetectorConstruction.hh"
#include "G4ThreeVector.hh"

class G4Material;
class G4VPhysicalVolume;
class DetectorConstruction;
class G4Material;

class DetectorConstructionSpec : public DetectorConstruction
{
public:
  DetectorConstructionSpec();
  ~DetectorConstructionSpec();

private:
  DetectorConstructionSpec( const DetectorConstructionSpec & );
  DetectorConstructionSpec &
  operator = ( const DetectorConstructionSpec & );

public:

private:
  virtual G4VPhysicalVolume *ConstructPayload( void );
  void MakeTrackers( G4VPhysicalVolume *pMother );
  void MakeCounters( G4VPhysicalVolume *pMother );
  void MakeSpecMagnetandTarget( G4VPhysicalVolume *pMother );
  void MakeVDetector( G4VPhysicalVolume *pMother );

  void SetRealMaterials( void );
  void PrintRealMaterialName( void ) const; 
  
};

#endif


















