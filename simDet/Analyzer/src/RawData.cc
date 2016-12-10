/*
  RawData.cc

  2016/4  K.Shirotori
*/

//Added Comments H. Asano
//role of this class
//
//1. decodes text files generated from GEANT4 simulation
//2. store Primary, T0 and SFT data into vectors
//The SFT data is stored in the vector of SFTRawHit class. 
//SFTRawHit is data storage object for SFT.

#include "RawData.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "Random/Randomize.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>
#include <string.h>

#include "DataType.hh"
#include "DetectorID.hh"
#include "PrimInfo.hh"
#include "HodoRawHit.hh"
#include "SFTRawHit.hh"
#include "TemplateLib.hh"

#include "ConfMan.hh"
#include "TrGeomMan.hh"

#define check1 0
#define check2 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
  PrimHC(0),
  T0RHC(0),
  SFTRawHitContainerArray() //vector of data object class for a single track in SFT
{}

RawData::~RawData()
{
  clearAll();
}

bool RawData::AddPrimInfo( PrimInfoContainer& cont,
			   double VertX, double VertY, double VertZ,
			   double BeamMom, double BeamU, double BeamV,
			   double AnaBeamMom, 
			   double Mass1, double Mom1,
			   double Theta1, double Phi1, double ThetaCM1, double PhiCM1,
			   double Mass2, double Mom2,
			   double Theta2, double Phi2, double ThetaCM2, double PhiCM2 )
{
  static const std::string funcname = "[RawData::AddPrimInfo]";
  
  PrimInfo *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    PrimInfo *q=cont[i];
    p=q; break;
  }
  if(!p){
    p = new PrimInfo();
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetVertX( VertX );
    p->SetVertY( VertY );
    p->SetVertZ( VertZ );
    p->SetBeamMom( BeamMom );
    p->SetBeamU( BeamU );
    p->SetBeamV( BeamV );
    p->SetAnaBeamMom( AnaBeamMom );
    p->SetMass1( Mass1 );
    p->SetMom1( Mom1 );
    p->SetTheta1( Theta1 );
    p->SetPhi1( Phi1 );
    p->SetThetaCM1( ThetaCM1 );
    p->SetPhiCM1( PhiCM1 );
    p->SetMass2( Mass2 );
    p->SetMom2( Mom2 );
    p->SetTheta2( Theta2 );
    p->SetPhi2( Phi2 );
    p->SetThetaCM2( ThetaCM2 );
    p->SetPhiCM2( PhiCM2 );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}


bool RawData::AddRawHit( SFTRawHitContainer& cont,
			 int Layer, int Ch, 
			 double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddRawHit]";
 
  SFTRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    SFTRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
	q->ChId()==Ch ){
      p=q; break;
    }
  }
  if(!p){
    p = new SFTRawHit( Layer, Ch );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetPosX( PosX );
    p->SetPosY( PosY );
    p->SetDL( DL );
    
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}


bool RawData::AddRawHit( SFTRawHitContainer& cont,
			 int Layer, int Ch, int hitID)
{
  static const std::string funcname = "[RawData::AddRawHit]";
 
  SFTRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    SFTRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
	q->ChId()==Ch ){
      p=q; break;
    }
  }
  if(!p){
    p = new SFTRawHit( Layer, Ch );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetHitID( hitID );
    
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddHodoRHit( HodoRHitContainer& cont,
			   int DetId, int Layer, int Seg,
			   double Time, double Edep, 
			   double Path, double Mom,
			   double PosX, double PosY, 
			   int Pid, double Beta )
{
  static const std::string funcname = "[RawData::AddHodoRHit]";

  HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->LayerId()==Layer &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit( DetId, Layer, Seg );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetTime( Time );
    p->SetEdep( Edep );
    p->SetPath( Path );
    p->SetMom( Mom );
    p->SetPosX( PosX );
    p->SetPosY( PosY );
    p->SetPid( Pid );
    p->SetBeta( Beta );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

void RawData::clearAll()
{
  std::for_each(PrimHC.begin(), PrimHC.end(), DeleteObject());
  PrimHC.clear();

  std::for_each(T0RHC.begin(), T0RHC.end(), DeleteObject());
  T0RHC.clear();

  for( int l=0; l<=PlMaxSFT; ++l){
    for_each( SFTRawHitContainerArray[l].begin(),  SFTRawHitContainerArray[l].end(), DeleteObject());
    SFTRawHitContainerArray[l].clear();
  }

  return;
}

//Decode data from text file 
//called event by event in EventBeamTracking::ProcessingNormal(std::ifstream &In)
bool RawData::DecodeSFTRawHits( std::ifstream &In )
{
  clearAll();
  ConfMan *confMan = ConfMan::GetConfManager();
  const TrGeomMan & geomMan=TrGeomMan::GetInstance();

#if check1
  std::cout << "***********************" << std::endl;
#endif

  int evNo=0;
  double type=0;
  In >> evNo;

#if check1
  std::cout << "ev = " << evNo << std::endl;
#endif
  
  int rawhitID = 0;
  while( type!=EventEnd ){
    In >> type;
    if( evNo==0 && type==0 ) return false;
#if check1
    std::cout << "type_0 = " << type << std::endl;
#endif


    ////////////////////////////////////////////////////////////////////////////////    
    ////PrimInfo
    ////////////////////////////////////////////////////////////////////////////////

    if( type==PrimInfoF ){
      while( type!=PrimInfoFEnd ){ 
	double x=0, y=0, z=0, pb=0, ub=0, vb=0, abmom=0;
	double m1=0, p1=0, theta1=0, phi1=0, thetacm1=0, phicm1=0;
	double m2=0, p2=0, theta2=0, phi2=0, thetacm2=0, phicm2=0;
	double pbeam=0; 
	
	In >> x;
	if( x==PrimInfoFEnd ) break;
	In >> y >> z 
	   >> pb >> ub >> vb >> abmom 
	   >> m1 >> p1 >> theta1 >> phi1 >> thetacm1 >> phicm1
	   >> m2 >> p2 >> theta2 >> phi2 >> thetacm2 >> phicm2;

	pbeam = abmom + CLHEP::RandGauss::shoot( 0.0, 0.0 );	

	AddPrimInfo( PrimHC, x, y, z, pb, ub, vb, pbeam, 
		     m1, p1, theta1, phi1, thetacm1, phicm1,
		     m2, p2, theta2, phi2, thetacm2, phicm2 );
	
#if check1
	std::cout << x  << " " << y  << " " << z  << " " 
		  << pb << " " << ub << " " << ub << " " << abmom << " " 
		  << m1 << " " << p1 << " " 
		  << theta1 << " " << phi1 << " " << thetacm1 << " " << phicm1 << " "
		  << m2 << " " << p2 << " " 
		  << theta2 << " " << phi2 << " " << thetacm2 << " " << phicm2 << " "
		  << std::endl;
#endif
      }
    }

#if check1
      std::cout << "type1_1 =" << type << std::endl;
#endif
    
    //dummy , when layer number starts from 1 in the previous code
    //AddTrRHit(SFTRawHitContainerArray[0], 0, 0, 0., 0., 0.);
    
    ////////////////////////////////////////////////////////////////////////////////
    ////Full tracking
    ////////////////////////////////////////////////////////////////////////////////
    
    if( type==FullTrackF ){ 
      while( type!=FullTrackFEnd ){ //loop over all tracks from here
	In >> type;
	
#if check1
	std::cout << "type1_2 =" << type << std::endl;
#endif
	  
	//Tracker
	if( type==FullTrackTF ){ 
	  int lnum=-1, ch=-1;
	  double x=0, y=0, dl=0;


	  while( lnum!=FullTrackTFEnd ){
	    In >> lnum;
	    if( lnum==FullTrackTFEnd ) break;
      
      //Type 0 simple detector
	    if( confMan->AnaMode()==0 ){
	      In >> x >> y;
	      
	      //SFT
        //PlMinSFT = 0 , PLMaxSFT=11 , PlOffsSFT=0; defined in detectorID.hh
	      if( lnum>=PlMinSFT+PlOffsSFT && lnum<=PlMaxSFT+PlOffsSFT ){
		double angle = geomMan.GetTiltAngle( lnum );
		double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetSFTResol() );
		AddRawHit(SFTRawHitContainerArray[lnum-PlOffsSFT], lnum, ch, x, y, dl);
	      }else{
          std::cout << "Invalid SFT layer number !!  " << lnum << std::endl;
        }
	    }	    
      
      //Type 1 ,2, 3 realistic detector
	    if( confMan->AnaMode()>=1 ){
	      In >> ch;
        if(ch>=20){std::cout << "Invalid ch number" << ch << std::endl;}
	      x=0.0; y=0.0; dl=0.0;//
	      
	      //SFT only stores layer number and segment id
	      if( lnum>=PlMinSFT+PlOffsSFT && lnum<=PlMaxSFT+PlOffsSFT ){
		      AddRawHit(SFTRawHitContainerArray[lnum-PlOffsSFT], lnum, ch, rawhitID);
          //std::cout << "lnum " << lnum << "ch " << ch << std::endl;
          //SFTRawHit *praw = new SFTRawHit(lnum,ch);
          //praw->SetHitID(rawhitID++);
          rawhitID++;
        }else{
          std::cerr << __FILE__ << "  " << __LINE__ 
          << " invalid SFT data format !! " << std::endl;
        }
	    }	    

#if check1
	    std::cout << lnum <<" " << ch << " " << x << " " << y << std::endl;
#endif
	  }
	}
	
#if check1
	std::cout << "type1_3 =" << type << std::endl;
#endif     
	
	//Counter
	if( type==FullTrackCF ){ 
	  int lnum=0;
	  int hits=0;
	  int layer=0, seg=0, pid=0;
	  double time=0, edep=0, path=0, mom=0, x=0, y=0, beta=0;
	  
	  while( lnum!=FullTrackCFEnd ){
	    In >> lnum;
	    if( lnum==FullTrackCFEnd ) break;
	    In >> hits;
	    
	    if( hits!=0 ){
	      In >> layer >> seg
		 >> time >> edep >> path >> mom 
		 >> x >> y >> pid >> beta;
	      
	      //T0
	      if( lnum == geomMan.GetDetectorId("T0") ){
		double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetT0Resol() );
		AddHodoRHit(T0RHC, lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
	      }
	      
#if check1	    
	      std::cout << lnum  << " " << hits << " " 
			<< layer << " " << seg  << " "
			<< time  << " " << edep << " " << path << " " << mom  << " "
			<< x     << " " << y    << " " << pid  << " " << beta << std::endl;
#endif
	    }
	  }
	}
      }
    }
    
#if check1
    std::cout << "type1_4 =" << type << std::endl;
#endif

  }//Event End

  return true;
}

const PrimInfoContainer& RawData::GetPrimHC() const
{
  return PrimHC;
}

const HodoRHitContainer& RawData::GetT0RHC() const
{
  return T0RHC;
}

const SFTRawHitContainer & RawData::GetSFTRawHitContainer( int layer ) const
{
  if( layer<0 || layer>PlMaxSFT ){
    std::cerr << __FILE__ << "  " << __LINE__ << " invalid layer number !! " << layer << std::endl;
    layer=-1;
  }
  return SFTRawHitContainerArray[layer];
}
