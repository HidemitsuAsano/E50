/*
  PhysicsList.cc

  2012/5  K.Shirotori
*/

#include "globals.hh"
#include "PhysicsList.hh"
#include "Transportation.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4KL3DecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
#include <iomanip>

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4FastSimulationManagerProcess.hh"

#include "ConfMan.hh"

PhysicsList::PhysicsList()
  : G4VUserPhysicsList()
{
  defaultCutValue = 0.1*cm;
}

PhysicsList::~PhysicsList()
{
}

void PhysicsList::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructHeavyIon();
}

void PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  // gamma
  G4Gamma::GammaDefinition();
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void PhysicsList::ConstructLeptons()
{
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void PhysicsList::ConstructMesons()
{
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();

  //Particles
  G4DecayTable *decayTable;
  G4VDecayChannel *mode;
  G4ParticleDefinition *particle;

  //K- no decay 
  particle 
    = new G4ParticleDefinition("kaonM1",    0.493677*GeV,  5.352e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,        -321,
			       true,         -1.0,          NULL,
			       false,       "kaon");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaonM1", 1.0, 2, "pi-", "pi0" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //K+ no decay 
  particle 
    = new G4ParticleDefinition("kaonP1",    0.493677*GeV,  5.352e-14*MeV,    +1.*eplus,
			       0,              -1,             0,
			       1,              +1,             0,
			       "meson",         0,             0,         321,
			       true,         -1.0,          NULL,
			       false,       "kaon");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "kaonP1", 1.0, 2, "pi+", "pi0" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //pi- no decay 
  particle 
    = new G4ParticleDefinition("pionM1",    0.1395700*GeV,  2.5452e-14*MeV,    -1.*eplus,
			       0,              -1,             0,
			       2,              -2,             -1,
			       "meson",         0,             0,         -211,
			       true,         -1.0,          NULL,
			       false,       "pi");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "pionM1", 1.0, 2, "mu-", "anti_nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //pi+ no decay 
  particle 
    = new G4ParticleDefinition("pionP1",    0.1395700*GeV,  2.5452e-14*MeV,    +1.*eplus,
			       0,              -1,             0,
			       2,              +2,             -1,
			       "meson",         0,             0,          211,
			       true,         -1.0,          NULL,
			       false,       "pi");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "pionP1", 1.0, 2, "mu+", "nu_mu" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  ///D mezom
  G4double Life_D0 = 0.0;
  G4double Life_D  = 0.0;
  // G4double litfe_D0 = 0.4101e-3*ns;
  // G4double litfe_D  = 1.040e-3*ns;
  //D0 -> K+ pi- mode (K+ and pi- -> no decay)
  particle 
    = new G4ParticleDefinition("DmesonZ1",    1.86484*GeV,  1.616e-9*MeV,  0.,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,         421,
			       false,         Life_D0,      NULL,
			       false,       "D");
  
  decayTable =  new G4DecayTable();
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ1", 1.0, 2, "kaon+", "pi-" );
  mode  = new G4PhaseSpaceDecayChannel( "DmesonZ1", 1.0, 2, "kaonP1", "pionM1" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //D0 -> K+ pi- pi+ pi- mode (K+, pi- and pi+ -> no decay)
  particle 
    = new G4ParticleDefinition("DmesonZ2",    1.86484*GeV,  1.616e-9*MeV,  0.,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,         421,
			       false,        Life_D0,       NULL,
			       false,       "D");
  
  decayTable =  new G4DecayTable();
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ2", 1.0, 4, "kaon+", "pi-", "pi+", "pi-" );
  mode  = new G4PhaseSpaceDecayChannel( "DmesonZ2", 1.0, 4, "kaonP1", "pionM1", "pionP1", "pionM1" );
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ2", 1.0, 3, "kaonP1", "pionM1", "pi0");  
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ2", 1.0, 3, "kaonP1", "mu-", "anti_nu_mu" );  
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ2", 1.0, 3, "kaonP1", "e-", "anti_nu_e" );  

  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //D0 -> K- pi+ pi- pi+ mode (K-, pi+ and pi- -> no decay)
  particle 
    = new G4ParticleDefinition("DmesonZ3",    1.86484*GeV,  1.616e-9*MeV,  0.,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,         421,
			       false,        Life_D0,       NULL,
			       false,       "D");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "DmesonZ3", 1.0, 4, "kaonM1", "pionP1", "pionM1", "pionP1" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //D- -> K+ pi- pi- mode (K+ and pi- -> no decay)
  particle 
    = new G4ParticleDefinition("DmesonM1",    1.86962*GeV,  6.37e-10*MeV,  -1.*eplus,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,          -411,
			       false,        Life_D,        NULL,
			       false,       "D");
  
  decayTable =  new G4DecayTable();
  mode  = new G4PhaseSpaceDecayChannel( "DmesonM1", 1.0, 3, "kaonP1", "pionM1", "pionM1" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );

  //D0 -> K- pi+ mode (K- and pi+ -> no decay)
  particle 
    = new G4ParticleDefinition("DmesonZ4",    1.86484*GeV,  1.616e-9*MeV,  0.,
			       0,              -1,             0,
			       1,              -1,             0,
			       "meson",         0,             0,         421,
			       false,        Life_D0,       NULL,
			       false,       "D");
  
  decayTable =  new G4DecayTable();
  //mode  = new G4PhaseSpaceDecayChannel( "DmesonZ4", 1.0, 2, "kaon-", "pi+" );
  mode  = new G4PhaseSpaceDecayChannel( "DmesonZ4", 1.0, 2, "kaonM1", "pionP1" );
  decayTable->Insert( mode );
  particle->SetDecayTable( decayTable );
}

void PhysicsList::ConstructBaryons()
{
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}

void PhysicsList::ConstructHeavyIon()
{
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void PhysicsList::ConstructProcess()
{
  ConfMan *confMan = ConfMan::GetConfManager();
  int flag = confMan->PhysFlag();

  AddTransportationSks();

  if( GetFPhysProcEM(flag) ){
    ConstructEM();
    if( GetFPhysProcHD(flag) ) ConstructHadronic();
  }
  if( GetFPhysProcDCY(flag)  ) ConstructDecay();

  // G4StepLimiter *stepLimiter = new G4StepLimiter();
  // G4UserSpecialCuts *userCuts = new G4UserSpecialCuts();
  // theParticleIterator->reset();
  // while ((*theParticleIterator)()){
  //   G4ParticleDefinition* particle = theParticleIterator->value();
  //   G4ProcessManager* pmanager = particle->GetProcessManager();
  //   G4String particleName = particle->GetParticleName();
  //   pmanager ->AddDiscreteProcess(stepLimiter);
  //   pmanager ->AddDiscreteProcess(userCuts);
  // }
}

///////Transportation
void PhysicsList::AddTransportationSks()
{
  Transportation* theTransportationProcess= new Transportation();

  // loop over all particles in G4ParticleTable
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if ( !particle->IsShortLived() ) {
      // Add transportation process for all particles other than  "shortlived"
      if ( pmanager == 0) {
        // Error !! no process manager
        G4Exception("PhysicsList::AddTransportation : no process manager!");
      } 
      else {
        // add transportation with ordering = ( -1, "first", "first" )
        pmanager->AddProcess(theTransportationProcess);
        pmanager->SetProcessOrderingToFirst(theTransportationProcess,
                                            idxAlongStep);
        pmanager ->SetProcessOrderingToFirst(theTransportationProcess,
                                             idxPostStep);
      }
    }
    else {
      // shortlived particle case
    }
  }
}

///////Cut
void PhysicsList::SetCuts()
{
  // Suppress error message int case e/gamma/proton do not exist
  G4int temp = GetVerboseLevel();
  // Retrive verbose level
  SetVerboseLevel(temp);

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");

  SetCutsWithDefault();  
  //  theParticleIterator->reset();
  //  while( (*theParticleIterator)() ){
  //    G4ParticleDefinition *particle=theParticleIterator->value();
  //    particle->SetApplyCutsFlag( true );
    //////////////////////////////////////////////////////////////////////
    //    G4cout << particle->GetParticleName() << " ==> ApplyCutFlag = "
    //     << particle->GetApplyCutsFlag() << G4endl;
    //////////////////////////////////////////////////////////////////////
  //  }
}

///////EM
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

//#include "G4MultipleScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4UserSpecialCuts.hh"

void PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if( particleName == "gamma" ){
      pManager->AddProcess( new G4UserSpecialCuts(), -1, -1, 1 );    
      pManager->AddDiscreteProcess( new G4PhotoElectricEffect() );
      pManager->AddDiscreteProcess( new G4ComptonScattering() );
      pManager->AddDiscreteProcess( new G4GammaConversion() );
    }
    else if( particleName == "e-" ){
      pManager->AddProcess( new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),     -1, -1, 4 );
    }
    else if( particleName == "e+" ){
      pManager->AddProcess( new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess( new G4eMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4eIonisation(),         -1,  2, 3 );
      pManager->AddProcess( new G4eBremsstrahlung(),     -1, -1, 4 );
      pManager->AddProcess( new G4eplusAnnihilation(),    1, -1, 5 );
    }
    else if( particleName == "mu+" || particleName == "mu-" ){
      pManager->AddProcess( new G4UserSpecialCuts(),      -1, -1, 1 );
      pManager->AddProcess( new G4MuMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess( new G4MuIonisation(),         -1,  2, 3 ); 
      pManager->AddProcess( new G4MuBremsstrahlung(),     -1, -1, 4 );
      pManager->AddProcess( new G4MuPairProduction(),     -1, -1, 5 );
    }
    else if( !(particle->IsShortLived()) && particle->GetPDGCharge()!=0 &&
	     !( particleName=="chargedgeantino" 
		|| particleName=="antichargedgeantino") ){
      pManager->AddProcess(new G4UserSpecialCuts(),     -1, -1, 1 );
      pManager->AddProcess(new G4hMultipleScattering(), -1,  1, 2 );
      pManager->AddProcess(new G4hIonisation(),         -1,  2, 3 );
    }
  }
}

///////Hadron
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4LElastic.hh"

//Pion
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"

//Kaon
#include "G4KaonPlusInelasticProcess.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4LEKaonZeroLInelastic.hh"

//Nucleon
#include "G4ProtonInelasticProcess.hh"
#include "G4LEProtonInelastic.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LENeutronInelastic.hh"

void PhysicsList:: ConstructHadronic()
{
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticProcess->RegisterMe(theElasticModel);
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    //Pion      
    if( particleName == "pi+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionPlusInelasticProcess* theInelasticProcess 
	= new G4PionPlusInelasticProcess( "inelastic" );
      G4LEPionPlusInelastic* theLEInelasticModel 
	= new G4LEPionPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "pi-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4PionMinusInelasticProcess* theInelasticProcess 
	= new G4PionMinusInelasticProcess( "inelastic" );
      G4LEPionMinusInelastic* theLEInelasticModel 
	= new G4LEPionMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Kaon
    else if( particleName == "kaon+" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonPlusInelasticProcess* theInelasticProcess 
	= new G4KaonPlusInelasticProcess( "inelastic" );
      G4LEKaonPlusInelastic* theLEInelasticModel 
	= new G4LEKaonPlusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon-" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonMinusInelasticProcess* theInelasticProcess 
	= new G4KaonMinusInelasticProcess( "inelastic" );
      G4LEKaonMinusInelastic* theLEInelasticModel 
	= new G4LEKaonMinusInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0S" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroSInelasticProcess* theInelasticProcess 
	= new G4KaonZeroSInelasticProcess( "inelastic" );
      G4LEKaonZeroSInelastic* theLEInelasticModel 
	= new G4LEKaonZeroSInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "kaon0L" ) {
      pmanager->AddProcess( theElasticProcess );
      G4KaonZeroLInelasticProcess* theInelasticProcess 
	= new G4KaonZeroLInelasticProcess( "inelastic" );
      G4LEKaonZeroLInelastic* theLEInelasticModel 
	= new G4LEKaonZeroLInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    //Nucleon
    else if( particleName == "proton" ) {
      pmanager->AddProcess( theElasticProcess );
      G4ProtonInelasticProcess* theInelasticProcess 
	= new G4ProtonInelasticProcess( "inelastic" );
      G4LEProtonInelastic* theLEInelasticModel 
	= new G4LEProtonInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if( particleName == "neutron" ) {
      pmanager->AddProcess( theElasticProcess );
      G4NeutronInelasticProcess* theInelasticProcess 
	= new G4NeutronInelasticProcess( "inelastic" );
      G4LENeutronInelastic* theLEInelasticModel 
	= new G4LENeutronInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

///////Decay
#include "G4Decay.hh"

void PhysicsList::ConstructDecay()
{
  G4Decay *theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pManager = particle->GetProcessManager();
    if( theDecayProcess->IsApplicable(*particle) ){
      pManager->AddProcess( theDecayProcess );
      pManager->SetProcessOrdering( theDecayProcess, idxPostStep );
      pManager->SetProcessOrdering( theDecayProcess, idxAtRest );
    }
  }
}
