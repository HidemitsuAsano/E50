/*
  RawData.cc

  2016/2  K.Shirotori
*/

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
#include "TrRawHit.hh"
#include "s_BeamRawHit.hh"
#include "s_ScatRawHit.hh"
#include "TemplateLib.hh"

#include "ConfMan.hh"
#include "TrGeomMan.hh"

#define check1 0
#define check2 0

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

RawData::RawData():
  PrimHC(0),
  T0RHC(0), TOFRHC(0), ITOFRHC(0), PADRHC(0),
  RICHRHC(0), PID1RHC(0), PID2RHC(0), MFRHC(0), 
  VDRHC(0),
  BSFTRHC(), SSFTRHC(), AFTRHC(),
  IT1RHC(), IT2RHC(),
  ST1RHC(), ST2RHC(),
  s_BeamRHC(), s_ScatRHC()
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

bool RawData::AddTrRHit( TrRHitContainer& cont,
			 int Layer, int Wire, 
			 double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddTrRHit]";
 
  TrRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    TrRawHit *q=cont[i];
    if( q->LayerId()==Layer &&
	q->WireId()==Wire ){
      p=q; break;
    }
  }
  if(!p){
    p = new TrRawHit( Layer, Wire );
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

//For simple tracking
bool RawData::AddsBGeneRHit( s_BeamRHitContainer& cont, 
			     int TrackID, int Type,
			     double Mom, 
			     double MomX, double MomY, double MomZ,  
			     int Pid, 
			     double PosX, double PosY )
{ 
  static const std::string funcname = "[RawData::AddsBGeneRHit]";

  s_BeamRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_BeamRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_BeamRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetMom( Mom );
    p->SetMomX( MomX );
    p->SetMomY( MomY );
    p->SetMomZ( MomZ );
    p->SetPid( Pid );
    p->SetVertX( PosX );
    p->SetVertY( PosY );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsBTrRHit( s_BeamRHitContainer& cont,
			   int TrackID, int Type,
			   int Layer, int Wire, 
			   double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddsBeamTrRHit]";
  
  s_BeamRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_BeamRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_BeamRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsTrRHit( Layer, Wire, PosX, PosY, DL );
  
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsBHodoRHit( s_BeamRHitContainer& cont,
			     int TrackID, int Type,
			     int DetId, int Layer, int Seg,
			     double Time, double Edep, 
			     double Path, double Mom,
			     double PosX, double PosY, 
			     int Pid, double Beta )
{
  static const std::string funcname = "[RawData::AddsBeamRHit]";

  s_BeamRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_BeamRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_BeamRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsHodoRHit( DetId, Layer, Seg, 
		     Time, Edep, Path, Mom,
		     PosX, PosY, Pid, Beta );
    
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsSGeneRHit( s_ScatRHitContainer& cont, 
			     int TrackID, int Type,
			     double Mom, 
			     double MomX, double MomY, double MomZ,  
			     int Pid, 
			     double PosX, double PosY )
{ 
  static const std::string funcname = "[RawData::AddsSGeneRHit]";
  
  s_ScatRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_ScatRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_ScatRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetMom( Mom );
    p->SetMomX( MomX );
    p->SetMomY( MomY );
    p->SetMomZ( MomZ );
    p->SetPid( Pid );
    p->SetVertX( PosX );
    p->SetVertY( PosY );

    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsSTrRHit( s_ScatRHitContainer& cont,
			   int TrackID, int Type,
			   int Layer, int Wire, 
			   double PosX, double PosY, double DL )
{
  static const std::string funcname = "[RawData::AddsScatTrRHit]";
  
  s_ScatRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_ScatRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_ScatRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsTrRHit( Layer, Wire, PosX, PosY, DL );
  
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

bool RawData::AddsSHodoRHit( s_ScatRHitContainer& cont,
			     int TrackID, int Type,
			     int DetId, int Layer, int Seg,
			     double Time, double Edep, 
			     double Path, double Mom,
			     double PosX, double PosY, 
			     int Pid, double Beta )
{
  static const std::string funcname = "[RawData::AddsBeamRHit]";

  s_ScatRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_ScatRawHit *q=cont[i];
    if( q->TrackId()==TrackID &&
	q->TrackType()==Type ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_ScatRawHit( TrackID, Type );
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetsHodoRHit( DetId, Layer, Seg, 
		     Time, Edep, Path, Mom,
		     PosX, PosY, Pid, Beta );
    
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

  std::for_each(TOFRHC.begin(), TOFRHC.end(), DeleteObject());
  TOFRHC.clear();

  std::for_each(ITOFRHC.begin(), ITOFRHC.end(), DeleteObject());
  ITOFRHC.clear();

  std::for_each(PADRHC.begin(), PADRHC.end(), DeleteObject());
  PADRHC.clear();

  std::for_each(RICHRHC.begin(), RICHRHC.end(), DeleteObject());
  RICHRHC.clear();

  std::for_each(PID1RHC.begin(), PID1RHC.end(), DeleteObject());
  PID1RHC.clear();

  std::for_each(PID2RHC.begin(), PID2RHC.end(), DeleteObject());
  PID2RHC.clear();

  std::for_each(MFRHC.begin(), MFRHC.end(), DeleteObject());
  MFRHC.clear();

  std::for_each(VDRHC.begin(), VDRHC.end(), DeleteObject());
  VDRHC.clear();

  for( int l=0; l<=PlMaxBSFT; ++l){
    for_each( BSFTRHC[l].begin(),  BSFTRHC[l].end(), DeleteObject());
    BSFTRHC[l].clear();
  }

  for( int l=0; l<=PlMaxSSFT; ++l){
    for_each( SSFTRHC[l].begin(),  SSFTRHC[l].end(), DeleteObject());
    SSFTRHC[l].clear();
  }

  for( int l=0; l<=(PlMaxAFTR+PlMaxAFTL); ++l){
    for_each( AFTRHC[l].begin(),  AFTRHC[l].end(), DeleteObject());
    AFTRHC[l].clear();
  }

  for( int l=0; l<=PlMaxIT1; ++l){
    for_each( IT1RHC[l].begin(),  IT1RHC[l].end(), DeleteObject());
    IT1RHC[l].clear();
  }

  for( int l=0; l<=(PlMaxIT2R+PlMaxIT2L); ++l){
    for_each( IT2RHC[l].begin(),  IT2RHC[l].end(), DeleteObject());
    IT2RHC[l].clear();
  }

  for( int l=0; l<=PlMaxST1; ++l){
    for_each( ST1RHC[l].begin(),  ST1RHC[l].end(), DeleteObject());
    ST1RHC[l].clear();
  }

  for( int l=0; l<=PlMaxST2; ++l){
    for_each( ST2RHC[l].begin(),  ST2RHC[l].end(), DeleteObject());
    ST2RHC[l].clear();
  }
    
  std::for_each(s_BeamRHC.begin(), s_BeamRHC.end(), DeleteObject());
  s_BeamRHC.clear();

  std::for_each(s_ScatRHC.begin(), s_ScatRHC.end(), DeleteObject());
  s_ScatRHC.clear();

  return;
}

bool RawData::DecodeRawHits( std::ifstream &In )
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
	double x, y, z, pb, ub, vb, abmom;
	double m1, p1, theta1, phi1, thetacm1, phicm1;
	double m2, p2, theta2, phi2, thetacm2, phicm2;
	double pbeam; 
	
	In >> x;
	if( x==PrimInfoFEnd ) break;
	In >> y >> z 
	   >> pb >> ub >> vb >> abmom 
	   >> m1 >> p1 >> theta1 >> phi1 >> thetacm1 >> phicm1
	   >> m2 >> p2 >> theta2 >> phi2 >> thetacm2 >> phicm2;

	pbeam = abmom + CLHEP::RandGauss::shoot( 0.0, confMan->GetBeamResol() );	

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

    AddTrRHit(BSFTRHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(SSFTRHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(AFTRHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(IT1RHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(IT2RHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(ST1RHC[0], 0, 0, 0., 0., 0.);
    AddTrRHit(ST2RHC[0], 0, 0, 0., 0., 0.);

    
    ////////////////////////////////////////////////////////////////////////////////
    ////Full tracking
    ////////////////////////////////////////////////////////////////////////////////
    
    if( type==FullTrackF ){ 
      while( type!=FullTrackFEnd ){ 
	In >> type;
	
#if check1
	std::cout << "type1_2 =" << type << std::endl;
#endif
	  
	//Tracker
	if( type==FullTrackTF ){ 
	  int lnum, wire=1;
	  double x, y, dl;


	  while( lnum!=FullTrackTFEnd ){
	    In >> lnum;
	    if( lnum==FullTrackTFEnd ) break;
	    In >> x >> y;
	    
	    //BSFT
	    if( lnum>=PlMinBSFT+PlOffsBSFT && lnum<=PlMaxBSFT+PlOffsBSFT ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetBSFTResol() );
	      AddTrRHit(BSFTRHC[lnum-PlOffsBSFT], lnum, wire, x, y, dl);
	    }
	    
	    //SSFT
	    if( lnum>=PlMinSSFT+PlOffsSSFT && lnum<=PlMaxSSFT+PlOffsSSFT ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetSSFTResol() );
	      AddTrRHit(SSFTRHC[lnum-PlOffsSSFT], lnum, wire, x, y, dl);
	    }
	    
	    //AFTR
	    if( lnum>=PlMinAFTR+PlOffsAFT && lnum<=PlMaxAFTR+PlOffsAFT ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetAFTResol() );
	      AddTrRHit(AFTRHC[lnum-PlOffsAFT], lnum, wire, x, y, dl);
	    }
	    
	    //AFTL
	    if( lnum>=PlMinAFTL+PlOffsAFT+2 && lnum<=PlMaxAFTL+PlOffsAFT+2 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetAFTResol() );
	      AddTrRHit(AFTRHC[lnum-PlOffsAFT], lnum, wire, x, y, dl);
	    }
	    
	    //IT1
	    if( lnum>=PlMinIT1+PlOffsIT1 && lnum<=PlMaxIT1+PlOffsIT1 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT1Resol() );
	      AddTrRHit(IT1RHC[lnum-PlOffsIT1], lnum, wire, x, y, dl);
	    }
	    
	    //IT2R
	    if( lnum>=PlMinIT2R+PlOffsIT2 && lnum<=PlMaxIT2R+PlOffsIT2 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT2Resol() );
	      AddTrRHit(IT2RHC[lnum-PlOffsIT2], lnum, wire, x, y, dl);
	    }
	    
	    //IT2L
	    if( lnum>=PlMinIT2L+PlOffsIT2+10 && lnum<=PlMaxIT2L+PlOffsIT2+10 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT2Resol() );
	      AddTrRHit(IT2RHC[lnum-PlOffsIT2-1], lnum, wire, x, y, dl);
	    }
	    
	    //ST1
	    if( lnum>=PlMinST1+PlOffsST1 && lnum<=PlMaxST1+PlOffsST1 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetST1Resol() );
	      AddTrRHit(ST1RHC[lnum-PlOffsST1], lnum, wire, x, y, dl);
	    }
	    
	    //ST2
	    if( lnum>=PlMinST2+PlOffsST2 && lnum<=PlMaxST2+PlOffsST2 ){
	      double angle = geomMan.GetTiltAngle( lnum );
	      double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
	      dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetST2Resol() );
	      AddTrRHit(ST2RHC[lnum-PlOffsST2], lnum, wire, x, y, dl);
	    }
	    
#if check1
	    std::cout << lnum <<" " << wire << " " << x << " " << y << std::endl;
#endif
	  }
	}
	
#if check1
	std::cout << "type1_3 =" << type << std::endl;
#endif     
	
	//Counter
	if( type==FullTrackCF ){ 
	  int lnum;
	  int hits;
	  int layer, seg, pid;
	  double time, edep, path, mom, x, y, beta;
	  
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
	      
	      //TOF
	      if( lnum == geomMan.GetDetectorId("TOF") ){
		double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetTofResol() );
		AddHodoRHit(TOFRHC, lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //ITOF
	      if( lnum == geomMan.GetDetectorId("ITOFR") ||
		  lnum == geomMan.GetDetectorId("ITOFL") ){
		double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetITofResol1() );
		AddHodoRHit(ITOFRHC, lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //ITOFS
	      if( lnum == geomMan.GetDetectorId("ITOFSR") ||
		  lnum == geomMan.GetDetectorId("ITOFSL") ){
		double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetITofResol2() );
		AddHodoRHit(ITOFRHC, lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //PAD
	      if( lnum == geomMan.GetDetectorId("PADU") ||
		  lnum == geomMan.GetDetectorId("PADD") ){
		double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetPADResol() );
		AddHodoRHit(PADRHC, lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //RICH
	      if( lnum == geomMan.GetDetectorId("RICH") ){
		AddHodoRHit(RICHRHC, lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //PID1 , AC
	      if( lnum == geomMan.GetDetectorId("PID1R") ||
		  lnum == geomMan.GetDetectorId("PID1L") ){
		AddHodoRHit(PID1RHC, lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //PID2
	      if( lnum == geomMan.GetDetectorId("PID2R") ||
		  lnum == geomMan.GetDetectorId("PID2L") ){
		AddHodoRHit(PID2RHC, lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //MF
	      if( lnum == geomMan.GetDetectorId("MFU") ||
		  lnum == geomMan.GetDetectorId("MFD") ){
		AddHodoRHit(MFRHC, lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
	      }
	      
	      //VD
	      if( lnum == geomMan.GetDetectorId("SpecVp1") ||
		  lnum == geomMan.GetDetectorId("SpecVp2") ||
		  lnum == geomMan.GetDetectorId("SpecVp3") ){
		AddHodoRHit(VDRHC, lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
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
    

    ////////////////////////////////////////////////////////////////////////////////
    ////Simple tracking
    ////////////////////////////////////////////////////////////////////////////////

    int id=0;    
    int typ=0;
    if( type==SimpleTrackF ){ 
      while( type!=SimpleTrackFEnd ){ 
	In >> type;
	
	////Beam
	if( type==SimpleTrackBF ){
	  In >> id >> typ;
#if check2
	  std::cout << "Beam TrackID = " << id << std::endl;
	  std::cout << "Track Type = " << typ << std::endl;
#endif 
	  while( type!=SimpleTrackBFEnd ){ 
	    In >> type;
	
#if check2
	    std::cout << "type2_1= " << type << std::endl;
#endif 

	    //Generated info
	    if( type==SimpleTrackGF ){ 
	      while( type!=SimpleTrackGFEnd ){ 
		double p, px, py, pz, x, y;
		int pid;
		
		In >> p;
		if( p==SimpleTrackGFEnd ) break;
		In >> px >> py >> pz >> pid >> x >> y;

		AddsBGeneRHit(s_BeamRHC, id, typ, p, px, py, pz, pid, x, y);
		
#if check2
		std::cout << p   << " " << px  << " " << py << " " << pz  << " " 
			  << pid << " " << x   << " " << y  << std::endl;
#endif      
	      }
	    }
	    
#if check2
	    std::cout << "type2_2= " << type << std::endl;
#endif 

	    //Counter(T0) info
	    if( type==SimpleTrackCF ){ 
	      while( type!=SimpleTrackCFEnd ){ 
		int lnum;
		int layer, seg, pid;
		double time, edep, path, mom, x, y, beta;
		
		In >> lnum;
		if( lnum==SimpleTrackCFEnd ) break;
		In >> layer >> seg
		   >> time >> edep >> path >> mom 
		   >> x >> y >> pid >> beta;

		//T0
		if( lnum == geomMan.GetDetectorId("T0") ){
		  double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetT0Resol() );
		  AddsBHodoRHit(s_BeamRHC, id, typ,
				lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
		}

#if check2  
		std::cout << lnum  << " " << layer << " " << seg  << " "
			  << time  << " " << edep << " " << path << " " << mom  << " "
			  << x     << " " << y    << " " << pid  << " " << beta << std::endl;
#endif
	      }
	    }

#if check2
	    std::cout << "type2_3= " << type << std::endl;
#endif 

	    //Tracker(BSFT) info
	    if( type==SimpleTrackTF ){ 
	      while( type!=SimpleTrackTFEnd ){ 
		int lnum, wire=1;
		double x, y, dl;
		In >> lnum;
		if( lnum==SimpleTrackTFEnd ) break;
		
		In >> x >> y;

		//BSFT
		if( lnum>=PlMinBSFT+PlOffsBSFT && lnum<=PlMaxBSFT+PlOffsBSFT ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetBSFTResol() );
		  AddsBTrRHit(s_BeamRHC, id, typ, lnum, wire, x, y, dl);
		}
	    
#if check2
		std::cout << lnum << " " << x << " " << y << std::endl;
#endif
	      }
	    }
	  } 
	}

#if check2
	std::cout << "type2_4= " << type << std::endl;
#endif 

	////Scattered particles
	if( type==SimpleTrackSF ){
	  In >> id >> typ;
#if check2
	  std::cout << "Scattered TrackID = " << id << std::endl;
	  std::cout << "Track Type = " << typ << std::endl;
#endif 
	  while( type!=SimpleTrackSFEnd ){ 
	    In >> type;

#if check2
	    std::cout << "type2_5= " << type << std::endl;
#endif 
	    
	    //Generated info
	    if( type==SimpleTrackGF ){ 
	      while( type!=SimpleTrackGFEnd ){ 
		double p, px, py, pz, x, y;
		int pid;
		
		In >> p;
		if( p==SimpleTrackGFEnd ) break;
		In >> px >> py >> pz >> pid >> x >> y;

		AddsSGeneRHit(s_ScatRHC, id, typ, p, px, py, pz, pid, x, y);

#if check2
		std::cout << p   << " " << px  << " " << py << " " << pz  << " " 
			  << pid << " " << x   << " " << y  << std::endl;
#endif      
	      }
	    }

#if check2
	    std::cout << "type2_6= " << type << std::endl;
#endif 

	    //Counter(TOF, ITOF, RICH, PID1, PID2, MF) info
	    if( type==SimpleTrackCF ){ 
	      while( type!=SimpleTrackCFEnd ){ 
		int lnum;
		int layer, seg, pid;
		double time, edep, path, mom, x, y, beta;

		In >> lnum;
		if( lnum==SimpleTrackCFEnd ) break;
		In >> layer >> seg
		   >> time >> edep >> path >> mom 
		   >> x >> y >> pid >> beta;

		//TOF
		if( lnum == geomMan.GetDetectorId("TOF") ){
		  double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetTofResol() );
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
		}
		
		//ITOF
		if( lnum == geomMan.GetDetectorId("ITOFR") ||
		    lnum == geomMan.GetDetectorId("ITOFL") ){
		  double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetITofResol1() );
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
		}
		
		//ITOFS
		if( lnum == geomMan.GetDetectorId("ITOFSR") ||
		    lnum == geomMan.GetDetectorId("ITOFSL") ){
		  double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetITofResol2() );
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
		}
		
		//PAD
		if( lnum == geomMan.GetDetectorId("PADU") ||
		    lnum == geomMan.GetDetectorId("PADD") ){
		  double T = time + CLHEP::RandGauss::shoot( 0.0, confMan->GetPADResol() );
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, T, edep, path, mom, x, y, pid, beta);
		}
		
		//RICH
		if( lnum == geomMan.GetDetectorId("RICH") ){
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
		}
		//std::cout << "RICH hit " << x << "  " << y << std::endl;
		//PID1
		if( lnum == geomMan.GetDetectorId("PID1R") ||
		    lnum == geomMan.GetDetectorId("PID1L") ){
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
		}
		
		//PID2
		if( lnum == geomMan.GetDetectorId("PID2R") ||
		    lnum == geomMan.GetDetectorId("PID2L") ){
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
		}
		
		//MF
		if( lnum == geomMan.GetDetectorId("MFU") ||
		    lnum == geomMan.GetDetectorId("MFD") ){
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
		}
		
		//VD
		if( lnum == geomMan.GetDetectorId("SpecVp1") ||
		    lnum == geomMan.GetDetectorId("SpecVp2") ||
		    lnum == geomMan.GetDetectorId("SpecVp3") ){
		  AddsSHodoRHit(s_ScatRHC, id, typ,
				lnum, layer, seg, time, edep, path, mom, x, y, pid, beta);
		}

#if check2  
		std::cout << lnum  << " " << layer << " " << seg  << " "
			  << time  << " " << edep << " " << path << " " << mom  << " "
			  << x     << " " << y    << " " << pid  << " " << beta << std::endl;
#endif
	      }
	    }

#if check2
	    std::cout << "type2_7= " << type << std::endl;
#endif 

	    //Tracker(SSFT, IT1, IT2, ST1, ST2) info
	    if( type==SimpleTrackTF ){ 
	      while( type!=SimpleTrackTFEnd ){ 
		int lnum, wire=1;
		double x, y, dl;
		In >> lnum;
		if( lnum==SimpleTrackTFEnd ) break;
		In >> x >> y;

		//SSFT
		if( lnum>=PlMinSSFT+PlOffsSSFT && lnum<=PlMaxSSFT+PlOffsSSFT ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetSSFTResol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//AFTR
		if( lnum>=PlMinAFTR+PlOffsAFT && lnum<=PlMaxAFTR+PlOffsAFT ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetAFTResol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//AFTL
		if( lnum>=PlMinAFTL+PlOffsAFT+2 && lnum<=PlMaxAFTL+PlOffsAFT+2 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetAFTResol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//IT1
		if( lnum>=PlMinIT1+PlOffsIT1 && lnum<=PlMaxIT1+PlOffsIT1 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT1Resol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//IT2R
		if( lnum>=PlMinIT2R+PlOffsIT2 && lnum<=PlMaxIT2R+PlOffsIT2 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT2Resol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//IT2L
		if( lnum>=PlMinIT2L+PlOffsIT2+10 && lnum<=PlMaxIT2L+PlOffsIT2+10 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetIT2Resol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//ST1
		if( lnum>=PlMinST1+PlOffsST1 && lnum<=PlMaxST1+PlOffsST1 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetST1Resol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}
		
		//ST2
		if( lnum>=PlMinST2+PlOffsST2 && lnum<=PlMaxST2+PlOffsST2 ){
		  double angle = geomMan.GetTiltAngle( lnum );
		  double l = x*cos(angle*Deg2Rad) + y*sin(angle*Deg2Rad);
		  dl = l + CLHEP::RandGauss::shoot( 0.0, confMan->GetST2Resol() );
		  AddsSTrRHit(s_ScatRHC, id, typ, lnum, wire, x, y, dl);
		}

#if check2
		std::cout << lnum << " " << x << " " << y << std::endl;
#endif
	      }
	    }
	  } 
	}

#if check2
	std::cout << "type2_8= " << type << std::endl;
#endif 
      }
    }
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

const HodoRHitContainer& RawData::GetTOFRHC() const
{
  return TOFRHC;
}

const HodoRHitContainer& RawData::GetITOFRHC() const
{
  return ITOFRHC;
}

const HodoRHitContainer& RawData::GetPADRHC() const
{
  return PADRHC;
}

const HodoRHitContainer& RawData::GetRICHRHC() const
{
  return RICHRHC;
}

const HodoRHitContainer& RawData::GetPID1RHC() const
{
  return PID1RHC;
}

const HodoRHitContainer& RawData::GetPID2RHC() const
{
  return PID2RHC;
}

const HodoRHitContainer& RawData::GetMFRHC() const
{
  return MFRHC;
}

const HodoRHitContainer& RawData::GetVDRHC() const
{
  return VDRHC;
}

const TrRHitContainer & RawData::GetBSFTRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxBSFT ) layer=0;
  return BSFTRHC[layer];
}

const TrRHitContainer & RawData::GetSSFTRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxSSFT ) layer=0;
  return SSFTRHC[layer];
}

const TrRHitContainer & RawData::GetAFTRHC( int layer ) const
{
  if( layer<0 || layer>(PlMaxAFTR+PlMaxAFTL) ) layer=0;
  return AFTRHC[layer];
}

const TrRHitContainer & RawData::GetIT1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxIT1 ) layer=0;
  return IT1RHC[layer];
}

const TrRHitContainer & RawData::GetIT2RHC( int layer ) const
{
  if( layer<0 || layer>(PlMaxIT2R+PlMaxIT2L) ) layer=0;
  return IT2RHC[layer];
}

const TrRHitContainer & RawData::GetST1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxST1 ) layer=0;
  return ST1RHC[layer];
}

const TrRHitContainer & RawData::GetST2RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxST2 ) layer=0;
  return ST2RHC[layer];
}

const s_BeamRHitContainer & RawData::GetsBeamRHC( void ) const
{
  return s_BeamRHC;
}

const s_ScatRHitContainer & RawData::GetsScatRHC( void ) const
{
  return s_ScatRHC;
}
