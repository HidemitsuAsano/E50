/*
  MainSpec.cc

  2012/5  K.Shirotori
*/

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "ConfMan.hh"

#ifdef G4VIS_USE
#include "VisManager.hh"
#endif

#include "RunAction.hh"
#include "DetectorConstructionSpec.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "EventActionSpec.hh"
#include "SteppingAction.hh"
#include "AnalysisSpec.hh"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"
#include "EvtGenBase/EvtStdlibRandomEngine.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenExternal/EvtExternalGenList.hh"
#include "EvtGenExternal/EvtExternalGenFactory.hh"

const int BufSize = 144;

FILE *fp_acc;
int   flagAcceptChange;

////////////////////////////////////////////////////
// Usage:
//
// Prog <macro> <confFile> <RooTFile> <OutputFile>
//
/////////////////////////////////////////////////////

int main( int argc, char **argv )
{
  //
  // Initialization on Random Number Generator
  // Initial Seed is got from Kernel Entropy Pool
  //
#if 1
  int initSeed=GetIntFromKernelEntropyPool()&0x7FFFFFFF;
  CLHEP::HepRandom::setTheSeed(initSeed);
  int startSeed=CLHEP::HepRandom::getTheSeed();
  G4cout << "*** Initial Seed = " << startSeed << G4endl;
  CLHEP::HepRandom::showEngineStatus();
#endif

  //Run manager
  G4RunManager *runMan = new G4RunManager();

  //Conf manager
  G4String confFile = "conf/sim.conf.default";
  G4String CommandName = argv[0];

  char buf[BufSize], buf1[BufSize];

  // Arguments check & Configuration file
  if( argc>=3 ){
    confFile = argv[2];
  }
  else {
    std::cout << "Analyzer Configuration File ? [" << confFile
              << "] ===> ";
    if( fgets(buf,BufSize,stdin) )
      if( sscanf(buf,"%s",buf1)==1 )
        confFile=buf1;
  }

  ConfMan *confManager = new ConfMan( confFile );
  confManager->Initialize();

  G4String histname;
  if( argc >= 4 ){
    histname = argv[3];
  }else{
    histname = "geant4_test.root";
  }
  
  //Detector Construction
  DetectorConstruction *detector = new DetectorConstructionSpec;
  runMan->SetUserInitialization( detector );

  //  G4cout<<"--------- Press return ---------"<< G4endl;
  //  getchar();

  //Physics list
  PhysicsList *physList = new PhysicsList();
  runMan->SetUserInitialization( physList );
  
#ifdef G4VIS_USE
  VisManager *visMan = new VisManager();
  visMan->Initialize();
#endif

  //Analysis
  Analysis *anaMan = new AnalysisSpec( histname );

  char *dataname;
  if(argc >= 5){
    dataname = argv[4];
    anaMan->SetDataFile( dataname );
  }

  G4cout<< "************" <<G4endl;
  EvtStdlibRandomEngine eng;
  EvtRandom::setRandomEngine((EvtRandomEngine*)&eng);
  
  EvtAbsRadCorr* radCorrEngine = 0;
  std::list<EvtDecayBase*> extraModels;
  
#ifdef EVTGEN_EXTERNAL
  EvtExternalGenList genList;
  radCorrEngine = genList.getPhotosModel();
  extraModels = genList.getListOfModels();
#endif

  EvtExternalGenFactory* externalGenerators = EvtExternalGenFactory::getInstance();
  std::string xmlDir("./xmldoc");

  bool convertPhysCode(true);
  externalGenerators->definePythiaGenerator(xmlDir, convertPhysCode);

  EvtGen *evtgen(0);
  ConfMan *confMan = ConfMan::GetConfManager();
  if( confMan->ExistEvtData1() && confMan->ExistEvtData1()){
    //Initialize the generator - read in the decay table and particle properties
    evtgen = new EvtGen( confMan->EvtGenDecayName().c_str(),
			 confMan->EvtGenPDLName().c_str(),
			 (EvtRandomEngine*)&eng,
			 radCorrEngine, &extraModels);
  }
  G4cout<< "************" <<G4endl;

  //Primary generator
  PrimaryGeneratorAction *priGen =
    new PrimaryGeneratorAction( detector, anaMan, evtgen );
  //Run action
  RunAction *runAction = new RunAction( anaMan );
  //Event action
  EventAction *eventAction = new EventActionSpec( anaMan );
  //Stepping action
  SteppingAction* stepAction = new SteppingAction( detector, eventAction );

  runMan->SetUserAction( priGen );
  runMan->SetUserAction( runAction );
  runMan->SetUserAction( eventAction );
  runMan->SetUserAction( stepAction );
  
  runMan->Initialize();

  G4UImanager *UI = G4UImanager::GetUIpointer();
  if( argc==1 ){
    // Interactive Mode
    G4UItcsh     *UIshell = new G4UItcsh();
    G4UIsession  *session = new G4UIterminal( UIshell );
    UI->ApplyCommand("/control/execute prerun.mac");
    G4cout << "Begin Session \n";
    session->SessionStart();
    G4cout << "End Session " << G4endl;
    delete session;
  }
  else {
    // Batch Mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];

    UI->ApplyCommand(command+fileName);
  }

  delete anaMan;

#ifdef G4VIS_USE
  delete visMan;
#endif
 
  delete runMan;

  return 0;
}
