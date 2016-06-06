/*
  PGReaction6.hh

  2014/3  K.Shirotori
*/

#ifndef PGReaction6_h 
#define PGReaction6_h  1

class PrimaryGeneratorAction;
class G4Event;
class ReadDSDecayData;
class EvtGen;

class PGReaction6
{
public:
  PGReaction6( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction6(){}

  //pi- p -> Y*0 K*0 for Lambda
  void HyperonRKStarZFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Y*0 K*0 for Sigma0
  void HyperonRKStarZFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Y*- K*+ for Sigma-
  void HyperonRKStarPFromPiMinusProton3( G4Event *anEvent, EvtGen *evtGenerator );
  
private:
  PrimaryGeneratorAction *pGen;
};

#endif
