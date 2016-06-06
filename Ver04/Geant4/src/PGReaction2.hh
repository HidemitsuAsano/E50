/*
  PGReaction2.hh

  2013/7  K.Shirotori
*/

#ifndef PGReaction2_h 
#define PGReaction2_h  1

class PrimaryGeneratorAction;
class G4Event;
class EvtGen;

class PGReaction2
{
public:
  PGReaction2( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction2(){}

  //pi- p -> Lambdac+ D*-: FS All 
  void YcDStarMesonFromPiMinusProtonWithDecay1( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Yc*+ D*-: FS All 
  void YcDStarMesonFromPiMinusProtonWithDecay2( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Lambdac+ D**: FS All 
  void YcDStarMesonFromPiMinusProtonWithDecay3( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Lambdac+ D*- pi0: FS All 
  void YcDStarMesonFromPiMinusProtonWithDecay4( G4Event *anEvent, EvtGen *evtGenerator );

private:
  PrimaryGeneratorAction *pGen;
};

#endif
