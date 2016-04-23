/*
  PGReaction1.hh

  2013/4  K.Shirotori
*/

#ifndef PGReaction1_h 
#define PGReaction1_h  1

class PrimaryGeneratorAction;
class G4Event;
class ReadDSDecayData;
class EvtGen;

class PGReaction1
{
public:
  PGReaction1( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction1(){}

  //pi- p -> Yc+ D*- : FS K+ pi- pi- 
  void YcDStarMesonFromPiMinusProton1( G4Event *anEvent );
  //pi- p -> Yc+ D*- : FS All state  
  void YcDStarMesonFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Yc  D** : FS All state  
  void YcDStarMesonFromPiMinusProton3( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi- 
  void YcDStarMesonFromPiMinusProton4( G4Event *anEvent );
  //pi- p -> Yc+ D*- : FS K+ pi- pi-, Yc+ decay pi+ pi-
  void YcDStarMesonFromPiMinusProton5( G4Event *anEvent );
  //pi- p -> Yc+ D*- : FS K+ pi- pi+ pi- pi-, Yc+ decay pi+ pi-
  void YcDStarMesonFromPiMinusProton6( G4Event *anEvent );

private:
  PrimaryGeneratorAction *pGen;
};

#endif
