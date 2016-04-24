/*
  GeneratorHelper.cc
  
  2016/4  K.Shirotori
*/

#include "GeneratorHelper.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "ConfMan.hh"
#include "G4LorentzVector.hh"

#include <cmath>
#include <iomanip>

G4ThreeVector UniformDirectionInUV( G4double u0, G4double v0,
				    G4double hu, G4double hv )
{
  G4double du=0.,dv=0.;
  if(hu!=0.0) du=(G4UniformRand()-0.5)*hu;
  if(hv!=0.0) dv=(G4UniformRand()-0.5)*hv;
  
  G4double u=u0+du, v=v0+dv;
  G4double ninv=1./sqrt(1.+u*u+v*v);

  return G4ThreeVector( u*ninv, v*ninv, ninv );
}

G4ThreeVector GaussDirectionInUV( G4double u0, G4double v0,
				G4double su, G4double sv )
{
  G4double du=0.,dv=0.;
  if(su!=0.0) du=G4RandGauss::shoot(0.0,su);
  if(sv!=0.0) dv=G4RandGauss::shoot(0.0,sv);

  G4double u=u0+du, v=v0+dv;
  G4double ninv=1./sqrt(1.+u*u+v*v);

  return G4ThreeVector( u*ninv, v*ninv, ninv );
}

G4ThreeVector UniformPosition( G4double hx, G4double hy, G4double hz )
{
  G4double x=0., y=0., z=0.;
  if(hx!=0.0) x+=(G4UniformRand()-0.5)*hx;
  if(hy!=0.0) y+=(G4UniformRand()-0.5)*hy;
  if(hz!=0.0) z+=(G4UniformRand()-0.5)*hz;

  return G4ThreeVector( x, y, z );
}

G4ThreeVector GaussPosition( G4double sx, G4double sy, G4double hz )
{
  G4double x=0., y=0., z=0.;
  if(sx!=0.0) x+=G4RandGauss::shoot(0.0,sx);;
  if(sy!=0.0) y+=G4RandGauss::shoot(0.0,sy);;
  if(hz!=0.0) z+=(G4UniformRand()-0.5)*hz;

  return G4ThreeVector( x, y, z );

}

const G4double X0def=0.0, dXdef=0.0;    // mm
const G4double Y0def=0.0, dYdef=0.0;    // mm
const G4double Z0def=0.0, dZdef=0.0;    // mm
const G4double U0def=0.0, dUdef=0.0;    // no dim.
const G4double V0def=0.0, dVdef=0.0;    // no dim.
const G4double P0def=0.7, dPdef=0.0;    // GeV/c

BeamParam::BeamParam( PrimaryGeneratorAction *gact ) 
  : pid(1), x0(X0def*mm), y0(Y0def*mm), z0(Z0def*mm), u0(U0def), v0(V0def),
    p0(P0def*GeV), dx(dXdef*mm), dy(dYdef*mm), dz(dZdef*mm),
    du(dUdef), dv(dVdef), dp(dPdef*GeV),
    gun(0)
{
  // set paramaters from Conf Manager
  ConfMan *confMan = ConfMan::GetConfManager();
  pid=confMan->GetPidBeam();
  x0=confMan->GetBeamX0()*mm; y0=confMan->GetBeamY0()*mm;
  z0=confMan->GetBeamZ0()*mm;
  u0=confMan->GetBeamU0(); v0=confMan->GetBeamV0();
  G4double val=confMan->GetBeamP0();
  if(val>=0.0) p0=val*GeV;
  val=confMan->GetBeamDX(); if(val>=0.0) dx=val*mm;
  val=confMan->GetBeamDY(); if(val>=0.0) dy=val*mm;
  val=confMan->GetBeamDZ(); if(val>=0.0) dz=val*mm;
  val=confMan->GetBeamDU(); if(val>=0.0) du=val;
  val=confMan->GetBeamDV(); if(val>=0.0) dv=val;
  val=confMan->GetBeamDP(); if(val>=0.0) dp=val*GeV;

  gun=gact->chooseGun(pid);
}

