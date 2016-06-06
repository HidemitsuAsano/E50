/*
  GeneratorHelper.cc
  
  2012/5  K.Shirotori
*/

#include "GeneratorHelper.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "ConfMan.hh"
#include "G4LorentzVector.hh"
#include "AngDisGenerator.hh"

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

G4ThreeVector UniformDirectionInThetaPhi( G4double cost1, G4double cost2 )
{
  AGUniform gen( cost1, cost2 );
  return gen.GenerateDirection();
}

G4double BreitWigner( G4double mean, G4double gamma )
{
  G4double r=G4UniformRand()-0.5;
  return mean+0.5*gamma*tan(acos(-1.)*r);
}


const G4double X0def=0.0, dXdef=0.0;    // mm
const G4double Y0def=0.0, dYdef=0.0;    // mm
const G4double Z0def=0.0, dZdef=0.0;    // mm
const G4double U0def=0.0, dUdef=0.0;    // no dim.
const G4double V0def=0.0, dVdef=0.0;    // no dim.
const G4double P0def=0.7, dPdef=0.0;    // GeV/c
const G4double Theta1=0.0, Theta2=20.0; // degree

BeamParam::BeamParam( DetectorConstruction *det,
		      PrimaryGeneratorAction *gact ) 
  : pid(4), x0(X0def*mm), y0(Y0def*mm), z0(Z0def*mm), u0(U0def), v0(V0def),
    p0(P0def*GeV), dx(dXdef*mm), dy(dYdef*mm), dz(dZdef*mm),
    du(dUdef), dv(dVdef), dp(dPdef*GeV), dztarg(dZdef*mm),
    cost1(cos(Theta1*degree)), cost2(cos(Theta2*degree)),
    gTargetPos(det->TargetPosition()), 
    TargetRotAngH(det->TargetAngle()), RecMass(1.0*GeV), gun(0)
{
  // set paramaters from Conf Manager
  ConfMan *confMan = ConfMan::GetConfManager();
  pid=confMan->GetPidBeam();
  x0=confMan->GetBeamX0()*mm; y0=confMan->GetBeamY0()*mm;
  z0=confMan->GetBeamZ0()*mm;
  u0=confMan->GetBeamU0(); v0=confMan->GetBeamV0();
  RecMass=confMan->GetRecoilMass()*GeV;
  RecWidth=confMan->GetRecoilWidth()*MeV;
  G4double val=confMan->GetBeamP0();
  if(val>=0.0) p0=val*GeV;
  val=confMan->GetBeamDX(); if(val>=0.0) dx=val*mm;
  val=confMan->GetBeamDY(); if(val>=0.0) dy=val*mm;
  val=confMan->GetBeamDZ(); if(val>=0.0) dz=val*mm;
  val=confMan->GetBeamDU(); if(val>=0.0) du=val;
  val=confMan->GetBeamDV(); if(val>=0.0) dv=val;
  val=confMan->GetBeamDP(); if(val>=0.0) dp=val*GeV;
  val=confMan->GetScatAngTheta1();
  if(val>=0.0 && val<=180.0) cost1=cos(val*degree);
  val=confMan->GetScatAngTheta2();
  if(val>=0.0 && val<=180.0) cost2=cos(val*degree);

  gun=gact->chooseGun(pid);

  if(dz!=0.0) dztarg=det->TargetLength();
}

bool Decay2Body( G4double Mini, G4double Mf1, G4double Mf2,
		 const G4ThreeVector & Pini,
		 G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		 G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM,
		 const AngDisGenerator & generator )
{
  if(Mini<Mf1+Mf2){
    // std::cerr << "Mini < Mf1+Mf2 Mini=" << Mini/GeV << "GeV/c2 "
    // 	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
    // 	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  
  G4ThreeVector beta( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  G4double Ecmf1=(Mini*Mini-Mf2*Mf2+Mf1*Mf1)/Mini*0.5;
  G4double Ecmf2=(Mini*Mini-Mf1*Mf1+Mf2*Mf2)/Mini*0.5;
  G4double Pcm=sqrt((Mini*Mini-(Mf1+Mf2)*(Mf1+Mf2))*
		  (Mini*Mini-(Mf1-Mf2)*(Mf1-Mf2)))/Mini*0.5;

  G4ThreeVector UnitDir=generator.GenerateDirection();

  UnitDir.rotateUz( Pini.unit() );
 
  //  std::cout << "UnitDir=" << UnitDir << " Pini.unit()=" 
  //	    << Pini.unit() << std::endl;


  G4ThreeVector Pcmf1 =  Pcm*UnitDir;
  G4ThreeVector Pcmf2 = -Pcm*UnitDir;

  G4LorentzVector LVf1( Pcmf1, Ecmf1 ), LVf2( Pcmf2, Ecmf2 );
  Pf1CM=LVf1.vect();
  Pf2CM=LVf2.vect();
 
  LVf1.boost(beta); Pf1=LVf1.vect(); 
  LVf2.boost(beta); Pf2=LVf2.vect();

  //  std::cout << "CosTCM=" << cost << " PhiCM=" << phi/degree 
  //	    << " degree" << std::endl;
  //  std::cout << "PCM=" << Pcm/GeV << " " << Pcmf1/GeV << " -->"
  //    	    << Pf1/GeV << " " << Pf1.mag()/GeV << std::endl;

  return true;
}

bool Scattering2Body( G4double Mi1, G4double Mi2, G4double Mf1, G4double Mf2,
		      const G4ThreeVector & Pini,
		      G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		      G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM,
		      const AngDisGenerator & generator )
{
  //  std::cout << "Pini=" << Pini/GeV << "GeV/c" << std::endl;
  
  G4double Mc=sqrt(Mi1*Mi1 + Mi2*Mi2 + 2.*Mi2*sqrt(Mi1*Mi1+Pini.mag2()) );
  return Decay2Body( Mc, Mf1, Mf2, Pini, Pf1, Pf2, Pf1CM, Pf2CM, generator );
}


bool Decay3BodyPhaseSpace( G4double Mini, 
			   G4double Mf1, G4double Mf2, G4double Mf3,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1,  
			   G4ThreeVector & Pf2,
			   G4ThreeVector & Pf3 )
{
  // from G4PhaseSpaceDecayChannell::ManyBodyDecayIt()

  G4int index, index2;
  G4int numberOfDaughters=3;

  G4double parentmass = Mini;
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  daughtermass[0] = Mf1; 
  daughtermass[1] = Mf2;
  daughtermass[2] = Mf3; 

  G4double sumofdaughtermass =Mf1+Mf2+Mf3;
  if(sumofdaughtermass>Mini){
    std::cerr << "Mini < Mf1+Mf2+Mf3 Mini=" << Mini/GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0;
  
  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index =1; index < numberOfDaughters -1; index++) rd[index] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index =1; index < numberOfDaughters -1; index++) {
      for(index2 = index+1; index2 < numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]){
          temp        = rd[index];
          rd[index]  = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index =0; index < numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
    }
    delete [] rd;
    
    //Calculate daughter momentum
    weight = 1.0;
    index =numberOfDaughters-1;
    daughtermomentum[index]= Pmx( sm[index-1],daughtermass[index-1], sm[index]);
    for(index =numberOfDaughters-2; index>=0; index--) {
      // calculate 
      daughtermomentum[index]= Pmx( sm[index],daughtermass[index], sm[index +1]);
      if(daughtermomentum[index] < 0.0) {
        // !!! illegal momentum !!!
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	
	G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
		    "can not create decay products", JustWarning,
		    "sum of daughter mass is larger than parent mass");
	
	return 0;   // Error detection
	
      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
      }
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ",
                  " Cannot decay ", JustWarning,
		  " Decay Kinematics cannot be calculated ");
      
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return 0;  // Error detection
    }
  } while ( weight > G4UniformRand());

  G4double costheta, sintheta, phi; 
  G4double beta;
  G4ThreeVector LPf[numberOfDaughters];
  G4LorentzVector LVf[numberOfDaughters];
  
  index = numberOfDaughters -2;
  //costheta = -1.0;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));

  LPf[index]=direction*daughtermomentum[index];
  LVf[index].setVect(LPf[index]);
  LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  LPf[index+1]=(-1.0)*direction*daughtermomentum[index];
  LVf[index+1].setVect(LPf[index+1]);
  LVf[index+1].setE(sqrt(LPf[index+1].mag2()+daughtermass[index+1]*daughtermass[index+1]));

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    //costheta = -1.0;
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (G4int index2 = index+1; index2<numberOfDaughters; index2++) {
      // make G4LorentzVector for secondaries
      G4LorentzVector p4;
      p4.setVect(LPf[index2]);
      p4.setE(sqrt(LPf[index2].mag2()+daughtermass[index2]*daughtermass[index2]));

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      LVf[index2].setVect(p4.vect());
      LVf[index2].setE(p4.e());
     }
    LPf[index]=(-1.0)*direction*daughtermomentum[index];
    LVf[index].setVect(LPf[index]);
    LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  }
  
  G4ThreeVector betai( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  LVf[0].boost(betai); Pf1=LVf[0].vect();
  LVf[1].boost(betai); Pf2=LVf[1].vect();
  LVf[2].boost(betai); Pf3=LVf[2].vect();
 
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return true;
}

// bool Decay3BodyPhaseSpace( G4double Mini, G4double Mf1, G4double Mf2, G4double Mf3,
// 			   const G4ThreeVector & Pini,
// 			   G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
// 			   G4ThreeVector & Pf3 )
// {
//   // from G4PhaseSpaceDecayChannell::ThreeBodyDecayIt()
//   // algorism of the code is originally written in GDECA3 of GEANT3

//   G4double sumOfDMass=Mf1+Mf2+Mf3;
//   if(sumOfDMass>Mini){
//     std::cerr << "Mini < Mf1+Mf2+Mf3 Mini=" << Mini/GeV << "GeV/c2 "
// 	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
// 	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
// 	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " << std::endl;
//     return false;
//   }

//   G4double Momf1, Momf2, Momf3;
//   G4double MomMax=0.0, MomSum=0.0;

//   do {
//     G4double rd1=G4UniformRand();
//     G4double rd2=G4UniformRand();
//     if( rd2>rd1 ){ G4double t=rd1; rd1=rd2; rd2=t; }
//     MomMax=0.0; MomSum=0.0;
//     // daughter 1
//     G4double ene=rd2*(Mini-sumOfDMass);
//     Momf1=sqrt(ene*ene+2.0*ene*Mf1);
//     if( Momf1>MomMax ) MomMax=Momf1;
//     MomSum += Momf1;
//     // daughter 2
//     ene=(1.-rd1)*(Mini-sumOfDMass);
//     Momf2=sqrt(ene*ene+2.0*ene*Mf2);
//     if( Momf2>MomMax ) MomMax=Momf2;
//     MomSum += Momf2;
//     // daughter 3
//     ene=(rd1-rd2)*(Mini-sumOfDMass);
//     Momf3=sqrt(ene*ene+2.0*ene*Mf3);
//     if( Momf2>MomMax ) MomMax=Momf2;
//     MomSum += Momf3;
//   } while( MomMax > MomSum-MomMax );

//   G4double cost=2.*G4UniformRand()-1.0;
//   G4double sint=sqrt(1.-cost*cost);
//   G4double phi=G4UniformRand()*360.*degree;
//   G4double sinp=sin(phi), cosp=cos(phi);
//   G4ThreeVector dir1 = G4ThreeVector( sint*cosp, sint*sinp, cost );
//   G4double costn=(Momf2*Momf2-Momf3*Momf3-Momf1*Momf1)/(2.0*Momf3*Momf1);
//   G4double sintn=sqrt(1.-cost*cost);
//   G4double phin=G4UniformRand()*360.*degree;
//   G4double sinpn=sin(phin), cospn=cos(phin);
//   G4ThreeVector dir3 =
//     G4ThreeVector( sintn*cospn*cost*cosp-sintn*sinpn*sinp+costn*sint*cosp,
// 		   sintn*cospn*cost*sinp+sintn*sinpn*cosp+costn*sint*sinp,
// 		   -sintn*cospn*sint+costn*cost );
//   G4ThreeVector LPf1=dir1*Momf1;
//   G4ThreeVector LPf3=dir3*(Momf3/dir3.mag());
//   G4ThreeVector LPf2=-(LPf1+LPf3);

//   G4ThreeVector beta( Pini/sqrt(Mini*Mini+Pini.mag2()) );

//   G4LorentzVector LVf1( sqrt(LPf1.mag2()+Mf1*Mf1), LPf1 ); 
//   G4LorentzVector LVf2( sqrt(LPf2.mag2()+Mf2*Mf2), LPf2 ); 
//   G4LorentzVector LVf3( sqrt(LPf3.mag2()+Mf3*Mf3), LPf3 ); 

//   LVf1.boost(beta); Pf1=LVf1.vect();
//   LVf2.boost(beta); Pf2=LVf2.vect();
//   LVf3.boost(beta); Pf3=LVf3.vect();
 
//   return true;
// }

bool Decay3BodyPhaseSpaceCM( G4double Mini, G4double Mf1, G4double Mf2, G4double Mf3,
			     const G4ThreeVector & Pini,
			     G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
			     G4ThreeVector & Pf3, 
			     G4ThreeVector & Pf1CM,  G4ThreeVector & Pf2CM,
			     G4ThreeVector & Pf3CM )
{
  // from G4PhaseSpaceDecayChannell::ThreeBodyDecayIt()
  // algorism of the code is originally written in GDECA3 of GEANT3 
  G4double sumOfDMass=Mf1+Mf2+Mf3;
  if(sumOfDMass>Mini){
    std::cerr << "Mini < Mf1+Mf2+Mf3 Mini=" << Mini/GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  G4double Momf1, Momf2, Momf3;
  G4double MomMax=0.0, MomSum=0.0;

  do {
    G4double rd1=G4UniformRand();
    G4double rd2=G4UniformRand();
    if( rd2>rd1 ){ G4double t=rd1; rd1=rd2; rd2=t; }
    MomMax=0.0; MomSum=0.0;
    // daughter 1
    G4double ene=rd2*(Mini-sumOfDMass);
    Momf1=sqrt(ene*ene+2.0*ene*Mf1);
    if( Momf1>MomMax ) MomMax=Momf1;
    MomSum += Momf1;
    // daughter 2
    ene=(1.-rd1)*(Mini-sumOfDMass);
    Momf2=sqrt(ene*ene+2.0*ene*Mf2);
    if( Momf2>MomMax ) MomMax=Momf2;
    MomSum += Momf2;
    // daughter 3
    ene=(rd1-rd2)*(Mini-sumOfDMass);
    Momf3=sqrt(ene*ene+2.0*ene*Mf3);
    if( Momf2>MomMax ) MomMax=Momf2;
    MomSum += Momf3;
  } while( MomMax > MomSum-MomMax );

  G4double cost=2.*G4UniformRand()-1.0;
  G4double sint=sqrt(1.-cost*cost);
  G4double phi=G4UniformRand()*360.*degree;
  G4double sinp=sin(phi), cosp=cos(phi);
  G4ThreeVector dir1 = G4ThreeVector( sint*cosp, sint*sinp, cost );
  G4double costn=(Momf2*Momf2-Momf3*Momf3-Momf1*Momf1)/(2.0*Momf3*Momf1);
  G4double sintn=sqrt(1.-cost*cost);
  G4double phin=G4UniformRand()*360.*degree;
  G4double sinpn=sin(phin), cospn=cos(phin);
  G4ThreeVector dir3 =
    G4ThreeVector( sintn*cospn*cost*cosp-sintn*sinpn*sinp+costn*sint*cosp,
		   sintn*cospn*cost*sinp+sintn*sinpn*cosp+costn*sint*sinp,
		   -sintn*cospn*sint+costn*cost );
  G4ThreeVector LPf1=dir1*Momf1;
  G4ThreeVector LPf3=dir3*(Momf3/dir3.mag());
  G4ThreeVector LPf2=-(LPf1+LPf3);

  G4ThreeVector beta( Pini/sqrt(Mini*Mini+Pini.mag2()) );

  G4LorentzVector LVf1( sqrt(LPf1.mag2()+Mf1*Mf1), LPf1 ); 
  G4LorentzVector LVf2( sqrt(LPf2.mag2()+Mf2*Mf2), LPf2 ); 
  G4LorentzVector LVf3( sqrt(LPf3.mag2()+Mf3*Mf3), LPf3 ); 
  Pf1CM=LVf1.vect();
  Pf2CM=LVf2.vect();
  Pf3CM=LVf3.vect();

  LVf1.boost(beta); Pf1=LVf1.vect();
  LVf2.boost(beta); Pf2=LVf2.vect();
  LVf3.boost(beta); Pf3=LVf3.vect();
 
  return true;
}

bool Scattering3Body( G4double Mi1, G4double Mi2, 
		      G4double Mf1, G4double Mf2, G4double Mf3,
		      const G4ThreeVector & Pini,
		      G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
		      G4ThreeVector & Pf3,
		      const AngDisGenerator & generator )
{
  //  std::cout << "Pini=" << Pini/GeV << "GeV/c" << std::endl;

  G4double Mc=sqrt(Mi1*Mi1 + Mi2*Mi2 + 2.*Mi2*sqrt(Mi1*Mi1+Pini.mag2()) );
  return Decay3BodyPhaseSpace( Mc, Mf1, Mf2, Mf3, Pini, Pf1, Pf2, Pf3 );
}

bool Decay4BodyPhaseSpace( G4double Mini, 
			   G4double Mf1, G4double Mf2, 
			   G4double Mf3, G4double Mf4,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
			   G4ThreeVector & Pf3,  G4ThreeVector & Pf4 )
{
  // from G4PhaseSpaceDecayChannell::ManyBodyDecayIt()

  G4int index, index2;
  G4int numberOfDaughters=4;

  G4double parentmass = Mini;
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  daughtermass[0] = Mf1; 
  daughtermass[1] = Mf2;
  daughtermass[2] = Mf3; 
  daughtermass[3] = Mf4; 

  G4double sumofdaughtermass =Mf1+Mf2+Mf3+Mf4;
  if(sumofdaughtermass>Mini){
    std::cerr << "Mini < Mf1+Mf2+Mf3+Mf4 Mini=" << Mini/GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " 
	      << "Mf4=" <<  Mf4/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0;
  
  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index =1; index < numberOfDaughters -1; index++) rd[index] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index =1; index < numberOfDaughters -1; index++) {
      for(index2 = index+1; index2 < numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]){
          temp        = rd[index];
          rd[index]  = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index =0; index < numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
    }
    delete [] rd;
    
    //Calculate daughter momentum
    weight = 1.0;
    index =numberOfDaughters-1;
    daughtermomentum[index]= Pmx( sm[index-1],daughtermass[index-1], sm[index]);
    for(index =numberOfDaughters-2; index>=0; index--) {
      // calculate 
      daughtermomentum[index]= Pmx( sm[index],daughtermass[index], sm[index +1]);
      if(daughtermomentum[index] < 0.0) {
        // !!! illegal momentum !!!
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	
	G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
		    "can not create decay products", JustWarning,
		    "sum of daughter mass is larger than parent mass");
	
	return 0;   // Error detection
	
      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
      }
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ",
                  " Cannot decay ", JustWarning,
		  " Decay Kinematics cannot be calculated ");
      
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return 0;  // Error detection
    }
  } while ( weight > G4UniformRand());

  G4double costheta, sintheta, phi; 
  G4double beta;
  G4ThreeVector LPf[numberOfDaughters];
  G4LorentzVector LVf[numberOfDaughters];
  
  index = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));

  LPf[index]=direction*daughtermomentum[index];
  LVf[index].setVect(LPf[index]);
  LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  LPf[index+1]=(-1.0)*direction*daughtermomentum[index];
  LVf[index+1].setVect(LPf[index+1]);
  LVf[index+1].setE(sqrt(LPf[index+1].mag2()+daughtermass[index+1]*daughtermass[index+1]));

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (G4int index2 = index+1; index2<numberOfDaughters; index2++) {
      // make G4LorentzVector for secondaries
      G4LorentzVector p4;
      p4.setVect(LPf[index2]);
      p4.setE(sqrt(LPf[index2].mag2()+daughtermass[index2]*daughtermass[index2]));

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      LVf[index2].setVect(p4.vect());
      LVf[index2].setE(p4.e());
     }
    LPf[index]=(-1.0)*direction*daughtermomentum[index];
    LVf[index].setVect(LPf[index]);
    LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  }
  
  G4ThreeVector betai( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  LVf[0].boost(betai); Pf1=LVf[0].vect();
  LVf[1].boost(betai); Pf2=LVf[1].vect();
  LVf[2].boost(betai); Pf3=LVf[2].vect();
  LVf[3].boost(betai); Pf4=LVf[3].vect();
 
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return true;
}

bool Decay5BodyPhaseSpace( G4double Mini, 
			   G4double Mf1, G4double Mf2, 
			   G4double Mf3, G4double Mf4,
			   G4double Mf5,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
			   G4ThreeVector & Pf3,  G4ThreeVector & Pf4,
			   G4ThreeVector & Pf5 )
{
  // from G4PhaseSpaceDecayChannell::ManyBodyDecayIt()

  G4int index, index2;
  G4int numberOfDaughters=5;

  G4double parentmass = Mini;
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  daughtermass[0] = Mf1; 
  daughtermass[1] = Mf2;
  daughtermass[2] = Mf3; 
  daughtermass[3] = Mf4; 
  daughtermass[4] = Mf5; 

  G4double sumofdaughtermass =Mf1+Mf2+Mf3+Mf4+Mf5;
  if(sumofdaughtermass>Mini){
    std::cerr << "Mini < Mf1+Mf2+Mf3+Mf4+Mf5 Mini=" << Mini/GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " 
	      << "Mf4=" <<  Mf4/GeV << "GeV/c2 " 
	      << "Mf5=" <<  Mf5/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0;
  
  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index =1; index < numberOfDaughters -1; index++) rd[index] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index =1; index < numberOfDaughters -1; index++) {
      for(index2 = index+1; index2 < numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]){
          temp        = rd[index];
          rd[index]  = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index =0; index < numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
    }
    delete [] rd;
    
    //Calculate daughter momentum
    weight = 1.0;
    index =numberOfDaughters-1;
    daughtermomentum[index]= Pmx( sm[index-1],daughtermass[index-1], sm[index]);
    for(index =numberOfDaughters-2; index>=0; index--) {
      // calculate 
      daughtermomentum[index]= Pmx( sm[index],daughtermass[index], sm[index +1]);
      if(daughtermomentum[index] < 0.0) {
        // !!! illegal momentum !!!
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	
	G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
		    "can not create decay products", JustWarning,
		    "sum of daughter mass is larger than parent mass");
	
	return 0;   // Error detection
	
      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
      }
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ",
                  " Cannot decay ", JustWarning,
		  " Decay Kinematics cannot be calculated ");
      
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return 0;  // Error detection
    }
  } while ( weight > G4UniformRand());

  G4double costheta, sintheta, phi; 
  G4double beta;
  G4ThreeVector LPf[numberOfDaughters];
  G4LorentzVector LVf[numberOfDaughters];
  
  index = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));

  LPf[index]=direction*daughtermomentum[index];
  LVf[index].setVect(LPf[index]);
  LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  LPf[index+1]=(-1.0)*direction*daughtermomentum[index];
  LVf[index+1].setVect(LPf[index+1]);
  LVf[index+1].setE(sqrt(LPf[index+1].mag2()+daughtermass[index+1]*daughtermass[index+1]));

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (G4int index2 = index+1; index2<numberOfDaughters; index2++) {
      // make G4LorentzVector for secondaries
      G4LorentzVector p4;
      p4.setVect(LPf[index2]);
      p4.setE(sqrt(LPf[index2].mag2()+daughtermass[index2]*daughtermass[index2]));

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      LVf[index2].setVect(p4.vect());
      LVf[index2].setE(p4.e());
     }
    LPf[index]=(-1.0)*direction*daughtermomentum[index];
    LVf[index].setVect(LPf[index]);
    LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  }
  
  G4ThreeVector betai( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  LVf[0].boost(betai); Pf1=LVf[0].vect();
  LVf[1].boost(betai); Pf2=LVf[1].vect();
  LVf[2].boost(betai); Pf3=LVf[2].vect();
  LVf[3].boost(betai); Pf4=LVf[3].vect();
  LVf[4].boost(betai); Pf5=LVf[4].vect();
 
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return true;
}

bool Decay6BodyPhaseSpace( G4double Mini, 
			   G4double Mf1, G4double Mf2, 
			   G4double Mf3, G4double Mf4,
			   G4double Mf5, G4double Mf6,
			   const G4ThreeVector & Pini,
			   G4ThreeVector & Pf1,  G4ThreeVector & Pf2,
			   G4ThreeVector & Pf3,  G4ThreeVector & Pf4,
			   G4ThreeVector & Pf5,  G4ThreeVector & Pf6 )
{
  // from G4PhaseSpaceDecayChannell::ManyBodyDecayIt()

  G4int index, index2;
  G4int numberOfDaughters=6;

  G4double parentmass = Mini;
  G4double *daughtermass = new G4double[numberOfDaughters]; 
  daughtermass[0] = Mf1; 
  daughtermass[1] = Mf2;
  daughtermass[2] = Mf3; 
  daughtermass[3] = Mf4; 
  daughtermass[4] = Mf5; 
  daughtermass[5] = Mf6; 

  G4double sumofdaughtermass =Mf1+Mf2+Mf3+Mf4+Mf5+Mf6;
  if(sumofdaughtermass>Mini){
    std::cerr << "Mini < Mf1+Mf2+Mf3+Mf4+Mf5+Mf6 Mini=" << Mini/GeV << "GeV/c2 "
	      << "Mf1=" <<  Mf1/GeV << "GeV/c2 "
	      << "Mf2=" <<  Mf2/GeV << "GeV/c2 "
	      << "Mf3=" <<  Mf3/GeV << "GeV/c2 " 
	      << "Mf4=" <<  Mf4/GeV << "GeV/c2 " 
	      << "Mf5=" <<  Mf5/GeV << "GeV/c2 " 
	      << "Mf6=" <<  Mf6/GeV << "GeV/c2 " << std::endl;
    return false;
  }
  
  //Calculate daughter momentum
  G4double *daughtermomentum = new G4double[numberOfDaughters];
  G4ThreeVector direction;  
  G4double *sm = new G4double[numberOfDaughters];
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0;
  
  do {
    //Generate rundom number in descending order 
    G4double temp;
    G4double *rd = new G4double[numberOfDaughters];
    rd[0] = 1.0;
    for(index =1; index < numberOfDaughters -1; index++) rd[index] = G4UniformRand(); 
    rd[ numberOfDaughters -1] = 0.0;
    for(index =1; index < numberOfDaughters -1; index++) {
      for(index2 = index+1; index2 < numberOfDaughters; index2++) {
        if (rd[index] < rd[index2]){
          temp        = rd[index];
          rd[index]  = rd[index2];
          rd[index2] = temp;
        }
      }
    }
    
    //calcurate virtual mass 
    tmas = parentmass -  sumofdaughtermass;
    temp =  sumofdaughtermass; 
    for(index =0; index < numberOfDaughters; index++) {
      sm[index] = rd[index]*tmas + temp;
      temp -= daughtermass[index];
    }
    delete [] rd;
    
    //Calculate daughter momentum
    weight = 1.0;
    index =numberOfDaughters-1;
    daughtermomentum[index]= Pmx( sm[index-1],daughtermass[index-1], sm[index]);
    for(index =numberOfDaughters-2; index>=0; index--) {
      // calculate 
      daughtermomentum[index]= Pmx( sm[index],daughtermass[index], sm[index +1]);
      if(daughtermomentum[index] < 0.0) {
        // !!! illegal momentum !!!
	delete [] sm;
	delete [] daughtermass;
	delete [] daughtermomentum;
	
	G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt",
		    "can not create decay products", JustWarning,
		    "sum of daughter mass is larger than parent mass");
	
	return 0;   // Error detection
	
      } else {
	// calculate weight of this events
        weight *=  daughtermomentum[index]/sm[index];
      }
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ >100) {
      G4Exception("G4PhaseSpaceDecayChannel::ManyBodyDecayIt: ",
                  " Cannot decay ", JustWarning,
		  " Decay Kinematics cannot be calculated ");
      
      delete [] sm;
      delete [] daughtermass;
      delete [] daughtermomentum;
      return 0;  // Error detection
    }
  } while ( weight > G4UniformRand());

  G4double costheta, sintheta, phi; 
  G4double beta;
  G4ThreeVector LPf[numberOfDaughters];
  G4LorentzVector LVf[numberOfDaughters];
  
  index = numberOfDaughters -2;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  direction.setZ(costheta);
  direction.setY(sintheta*std::sin(phi));
  direction.setX(sintheta*std::cos(phi));

  LPf[index]=direction*daughtermomentum[index];
  LVf[index].setVect(LPf[index]);
  LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  LPf[index+1]=(-1.0)*direction*daughtermomentum[index];
  LVf[index+1].setVect(LPf[index+1]);
  LVf[index+1].setE(sqrt(LPf[index+1].mag2()+daughtermass[index+1]*daughtermass[index+1]));

  for (index = numberOfDaughters -3;  index >= 0; index--) {
    //calculate momentum direction
    costheta = 2.*G4UniformRand()-1.0;
    sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
    phi  = twopi*G4UniformRand()*rad;
    direction.setZ(costheta);
    direction.setY(sintheta*std::sin(phi));
    direction.setX(sintheta*std::cos(phi));

    beta = daughtermomentum[index];
    beta /= std::sqrt( daughtermomentum[index]*daughtermomentum[index] + sm[index+1]*sm[index+1] );
    for (G4int index2 = index+1; index2<numberOfDaughters; index2++) {
      // make G4LorentzVector for secondaries
      G4LorentzVector p4;
      p4.setVect(LPf[index2]);
      p4.setE(sqrt(LPf[index2].mag2()+daughtermass[index2]*daughtermass[index2]));

      // boost secondaries to  new frame 
      p4.boost( direction.x()*beta, direction.y()*beta, direction.z()*beta);

      // change energy/momentum
      LVf[index2].setVect(p4.vect());
      LVf[index2].setE(p4.e());
     }
    LPf[index]=(-1.0)*direction*daughtermomentum[index];
    LVf[index].setVect(LPf[index]);
    LVf[index].setE(sqrt(LPf[index].mag2()+daughtermass[index]*daughtermass[index]));
  }
  
  G4ThreeVector betai( Pini/sqrt(Mini*Mini+Pini.mag2()) );
  
  LVf[0].boost(betai); Pf1=LVf[0].vect();
  LVf[1].boost(betai); Pf2=LVf[1].vect();
  LVf[2].boost(betai); Pf3=LVf[2].vect();
  LVf[3].boost(betai); Pf4=LVf[3].vect();
  LVf[4].boost(betai); Pf5=LVf[4].vect();
  LVf[5].boost(betai); Pf6=LVf[5].vect();
 
  delete [] daughtermomentum;
  delete [] daughtermass;
  delete [] sm;
  
  return true;
}

G4double Pmx(G4double e, G4double p1, G4double p2)
{
   // calcurate momentum of daughter particles in two-body decay
   G4double ppp = (e+p1+p2)*(e+p1-p2)*(e-p1+p2)*(e-p1-p2)/(4.0*e*e);
   if (ppp>0) return std::sqrt(ppp);
   else       return -1.;
}
