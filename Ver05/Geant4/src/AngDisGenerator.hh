/*
  AngDisGenerator.hh

  2012/5  K.Shirotori
*/

#ifndef AngDisGenerarot_h
#define AngDisGenerarot_h 1

#include "G4ThreeVector.hh"

class AngDisGenerator
{
public:
  AngDisGenerator( G4double cost1=1.0, G4double cost2=-1.0 );
  virtual ~AngDisGenerator() {};

  virtual G4ThreeVector GenerateDirection( void ) const = 0;

protected:
  G4double cost1_, cost2_;

};

class AGSWave : public AngDisGenerator
{
  // D(x)=1/2 (Flat)
public:
  AGSWave( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGSWave() {};

  G4ThreeVector GenerateDirection( void ) const;
};

typedef AGSWave AGUniform;


class AGPWaveFP : public AngDisGenerator
{
  // D(x)=1/2*(1+x) (Forward p-wave)
public:
  AGPWaveFP( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGPWaveFP() {};

  G4ThreeVector GenerateDirection( void ) const;
};

class AGPWaveBP : public AngDisGenerator
{
  // D(x)=1/2*(1-x) (backward p-wave)
public:
  AGPWaveBP( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGPWaveBP() {};

  G4ThreeVector GenerateDirection( void ) const;
};

class AGDWave1 : public AngDisGenerator
{
  // D(x)= 1/2*(1+x*x) (Forward d-wave)
public:
  AGDWave1( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGDWave1() {};

  G4ThreeVector GenerateDirection( void ) const;
private:
  inline G4double Dfunc( G4double x ) const { return 0.5*(1.+x*x); }
};

class AGDWave2 : public AngDisGenerator
{
  // D(x)= 1/2*(1-x*x) (backward d-wave)
public:
  AGDWave2( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGDWave2() {};

  G4ThreeVector GenerateDirection( void ) const;
private:
  inline G4double Dfunc( G4double x ) const { return 1/2*(1.-x*x); }
};

class AGForward1 : public AngDisGenerator
{
  // D(x)= exp(a*(1-x))
public:
  AGForward1( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGForward1() {};

  G4ThreeVector GenerateDirection( void ) const;
private:
  //Fixed coefficient for t=-25 GeV^2 case
  //Not correctly calculating kinematics
  inline G4double Dfunc( G4double x ) const { return exp(-20.*(1.-x)); }
  //inline G4double Dfunc( G4double x ) const { return exp(-18.*(1.-x)); }
};

class AGForward2 : public AngDisGenerator
{
  // D(x)= exp(a*(1-x))
public:
  AGForward2( G4double cost1=1.0, G4double cost2=-1.0 );
  ~AGForward2() {};

  G4ThreeVector GenerateDirection( void ) const;
private:
  inline G4double Dfunc( G4double x ) const { return exp(-2.*(1.-x)); }
};

#endif
