//r_SFT_OpticalPhysics 
//
// May 10th 2016  
// Hidemitsu Asano

#ifndef r_SFT_OpticalPhysics_h
#define r_SFT_OpticalPhysics_h 1

#include "globals.hh"

#include "G4OpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VPhysicsConstructor.hh"

class r_SFT_OpticalPhysics : public G4VPhysicsConstructor
{
  public:

    r_SFT_OpticalPhysics(G4bool toggle=true);
    virtual ~r_SFT_OpticalPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    G4OpWLS* GetWLSProcess() {return theWLSProcess;}
    G4Cerenkov* GetCerenkovProcess() {return theCerenkovProcess;}
    G4Scintillation* GetScintillationProcess() {return theScintProcess;}
    G4OpAbsorption* GetAbsorptionProcess() {return theAbsorptionProcess;}
    G4OpRayleigh* GetRayleighScatteringProcess() {return theRayleighScattering;}
    G4OpBoundaryProcess* GetBoundaryProcess() { return theBoundaryProcess;}

    void SetNbOfPhotonsCerenkov(G4int);

private:

    G4OpWLS*             theWLSProcess;
    G4Cerenkov*          theCerenkovProcess;
    G4Scintillation*     theScintProcess;
    G4OpAbsorption*      theAbsorptionProcess;
    G4OpRayleigh*        theRayleighScattering;
    G4OpBoundaryProcess* theBoundaryProcess;
 
    G4bool AbsorptionOn;

};
#endif
