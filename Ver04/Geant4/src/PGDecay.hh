/*
  PGDecay.hh

  2012/5  K.Shirotori
*/

#ifndef PGDecay_h 
#define PGDecay_h 1

class PrimaryGeneratorAction;
class G4Event;

class PGDecay
{
public:
  PGDecay( PrimaryGeneratorAction *PGAction )
    : pGen(PGAction)
  {}
  ~PGDecay() {}

  void DMesonMinus3BodyDecayKPPiMPiM( G4Event *anEvent,
				      bool generatePiP = true );

  void DMesonZero2BodyDecayKPPiM( G4Event *anEvent,
				  bool generatePiP = true );

  void KMinus2BodyDecayMuMNu( G4Event *anEvent,
			      bool generateNeutrino = false );

  void KMinus2BodyDecayPiMPi0( G4Event *anEvent,
			       bool generatePi0 = true ); 

  void KMinus3BodyDecayPiMPiMPiP( G4Event *anEvent,
				  bool generatePiP = true );

  void KMinus3BodyDecayEMPi0Nu( G4Event *anEvent,
				bool generatePi0 = true,
				bool genetateNeutrino = false );

  void KMinus3BodyDecayMuMPi0Nu( G4Event *anEvent,
				 bool generatePi0 = true,
				 bool genetateNeutrino = false );

  void KMinus3BodyDecayPiMPi0Pi0( G4Event *anEvent,
				  bool generatePi0 = true );


private:
  PrimaryGeneratorAction *pGen;
};

#endif 
