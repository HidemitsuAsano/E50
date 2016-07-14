/*
  PhysicsList.hh

  2012/5  K.Shirotori
*/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class PhysicsList: public G4VUserPhysicsList
{
public:
  PhysicsList();
  ~PhysicsList();
  
protected:
  // Construct particle and physics
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructHeavyIon();
  
private:
  // these methods Construct physics processes and register them
  void AddTransportationSks();
  void ConstructDecay();
  void ConstructEM();
  void ConstructHadronic();
};

#endif

 
