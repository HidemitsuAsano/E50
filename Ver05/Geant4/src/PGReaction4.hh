/*
  PGReaction4.hh

  2013/4  K.Shirotori
*/

#ifndef PGReaction4_h 
#define PGReaction4_h  1

class PrimaryGeneratorAction;
class G4Event;
class ReadJAMData;

class PGReaction4
{
public:
  PGReaction4( PrimaryGeneratorAction * PGAction )
    : pGen(PGAction) 
  {}
  ~PGReaction4(){}

  //JAM BG
  void BackgroundFromPiMinusProtonJAM1( G4Event *anEvent, ReadJAMData *readJAMData );

private:
  PrimaryGeneratorAction *pGen;
};

#endif
