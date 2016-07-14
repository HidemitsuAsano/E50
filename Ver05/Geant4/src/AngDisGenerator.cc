/*
  AngDisGenerator.cc

  2012/5  K.Shirotori
*/

#include "AngDisGenerator.hh"
#include "Randomize.hh"

#include <cmath>

G4double URand(){ return  G4UniformRand(); }


AngDisGenerator::AngDisGenerator( G4double cost1, G4double cost2 )
  : cost1_(cost1), cost2_(cost2)
{}

AGSWave::AGSWave( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGSWave::GenerateDirection( void ) const
{
  G4double cost=cost1_+URand()*(cost2_-cost1_);
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGPWaveFP::AGPWaveFP( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGPWaveFP::GenerateDirection( void ) const
{
  G4double a=2./(cost2_-cost1_)/(2.+cost2_+cost1_);
  G4double cost=-1.+sqrt((cost1_+1.)*(cost1_+1.)+2.*URand()/a);
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}


AGPWaveBP::AGPWaveBP( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGPWaveBP::GenerateDirection( void ) const
{
  G4double a=2./(cost2_-cost1_)/(2.-cost2_-cost1_);
  G4double cost=1.-sqrt((cost1_-1.)*(cost1_-1.)-2.*URand()/a);
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGDWave1::AGDWave1( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGDWave1::GenerateDirection( void ) const
{
  G4double cost=(cost2_-cost1_)*URand()+cost1_;
  G4double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

AGDWave2::AGDWave2( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGDWave2::GenerateDirection( void ) const
{
  G4double cost=(cost2_-cost1_)*URand()+cost1_;
  G4double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

//For exp(2t)
AGForward1::AGForward1( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGForward1::GenerateDirection( void ) const
{
  G4double cost=(cost2_-cost1_)*URand()+cost1_;
  G4double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}

//For exp(2t)
AGForward2::AGForward2( G4double cost1, G4double cost2 )
  : AngDisGenerator( cost1, cost2 )
{}

G4ThreeVector AGForward2::GenerateDirection( void ) const
{
  G4double cost=(cost2_-cost1_)*URand()+cost1_;
  G4double p=URand();
  while ( p>Dfunc(cost) ){
    cost=(cost2_-cost1_)*URand()+cost1_;
    p=URand();
  }
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=URand()*acos(-1.)*2.;
  G4double cosp=cos(phi), sinp=sin(phi);

  return G4ThreeVector( sint*cosp, sint*sinp, cost );
}
