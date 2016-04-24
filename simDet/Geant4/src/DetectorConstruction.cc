/*
  DetectorConstruction.cc

  2016/4  K.Shirotori
*/

#include "DetectorConstruction.hh"
#include "MaterialList.hh" 
#include "ConfMan.hh"
#include "DCGeomMan.hh"

#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

DetectorConstruction::DetectorConstruction()
  : mList_(0)
{
}

DetectorConstruction::~DetectorConstruction()
{
  delete mList_;
}
MaterialList *DetectorConstruction::DefineMaterials()
{
  // not safe for exception 
  if(mList_) delete mList_;
  return new MaterialList();
}

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  mList_ = DefineMaterials();

  ConfMan *confMan = ConfMan::GetConfManager();

  G4VPhysicalVolume *world = ConstructPayload();

  return world;
}

