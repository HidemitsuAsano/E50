/*
  PGReaction3.hh

  2013/4  K.Shirotori
*/

#ifndef PGReaction3_h 
#define PGReaction3_h  1

class PrimaryGeneratorAction;
class G4Event;
class EvtGen;

class PGReaction3
{
public:
  PGReaction3( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction3(){}

  //pi- p -> D*- D*+ n : FS All state & All state
  void DStarMesonPairFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> M(ccbar) n : FS All state
  void DStarMesonPairFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator );

private:
  PrimaryGeneratorAction *pGen;
};

#endif
