/*
  PGReaction5.hh

  2013/12  K.Shirotori
*/

#ifndef PGReaction5_h 
#define PGReaction5_h  1

class PrimaryGeneratorAction;
class G4Event;
class ReadDSDecayData;
class EvtGen;

class PGReaction5
{
public:
  PGReaction5( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction5(){}

  //pi- p -> Xi*- K*0 : FS K+ K+ pi- 
  void XiKStarMesonFromPiMinusProton1( G4Event *anEvent, EvtGen *evtGenerator );
  //pi- p -> Xi*- K*0 from Y* production : FS K+ K+ pi- 
  void XiKStarMesonFromPiMinusProton2( G4Event *anEvent, EvtGen *evtGenerator );

  //K- p -> Xi*- K+ 
  void XiKStarMesonFromKMinusProton1( G4Event *anEvent, EvtGen *evtGenerator );
  //K- p -> Xi*- K*+ 
  void XiKStarMesonFromKMinusProton2( G4Event *anEvent, EvtGen *evtGenerator );
  //K- p -> Xi*0 K*0 
  void XiKStarMesonFromKMinusProton3( G4Event *anEvent, EvtGen *evtGenerator );

private:
  PrimaryGeneratorAction *pGen;
};

#endif
