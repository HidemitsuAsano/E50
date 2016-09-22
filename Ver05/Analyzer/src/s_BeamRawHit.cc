/*
  s_BeamRawHit.cc

  2016/2
*/

#include "s_BeamRawHit.hh"
#include "s_TrRawHit.hh"
#include "s_HodoRawHit.hh"
#include "ConfMan.hh"
#include "TrGeomMan.hh"
#include "DetectorID.hh"
#include "TemplateLib.hh"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <cstdlib>

s_BeamRawHit::~s_BeamRawHit()
{
  clearRegisteredHits();
}

void s_BeamRawHit::clearRegisteredHits( void )
{
  static const std::string funcname="[s_BeamRawHit::clearRegisteredHits]";

  std::for_each(sT0RHC.begin(), sT0RHC.end(), DeleteObject());
  sT0RHC.clear();

  for( int l=0; l<=PlMaxBSFT; ++l){
    for_each( sBSFTRHC[l].begin(),  sBSFTRHC[l].end(), DeleteObject());
    sBSFTRHC[l].clear();
  }
}

bool s_BeamRawHit::SetsTrRHit( int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_BeamRawHit::SetTrRHit]";

  //BSFT
  if( Layer>=PlMinBSFT+PlOffsBSFT && Layer<=PlMaxBSFT+PlOffsBSFT ){
    AddsTrRHit(sBSFTRHC[Layer-PlOffsBSFT], Layer, Wire, PosX, PosY, DL);
  }

  return true;
}

bool s_BeamRawHit::SetsHodoRHit( int DetId, int Layer, int Seg,
				 double Time, double Edep, 
				 double Path, double Mom,
				 double PosX, double PosY, 
				 int Pid, double Beta )
{
  static const std::string funcname="[s_BeamRawHit::SetHodoRHit]";
  
  const TrGeomMan & geomMan=TrGeomMan::GetInstance();
  
  //T0
  if( DetId == geomMan.GetDetectorId("T0") ){
    AddsHodoRHit(sT0RHC,
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }

  return true;
}

bool s_BeamRawHit::AddsTrRHit( s_TrRHitContainer &cont,
			       int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_BeamRawHit::AddTrRHit]";
  
  s_TrRawHit *p=0;
  int nh=cont.size();

  for( int i=0; i<nh; ++i ){
    s_TrRawHit *q=cont[i];
    if( q->LayerId()==Layer && 
	q->WireId()==Wire ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_TrRawHit( Layer, Wire );
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

bool s_BeamRawHit::AddsHodoRHit( s_HodoRHitContainer &cont,
				 int DetId, int Layer, int Seg,
				 double Time, double Edep, 
				 double Path, double Mom,
				 double PosX, double PosY, 
				 int Pid, double Beta )
{
  static const std::string funcname="[s_BeamRawHit::AddHodoRHit]";
  
  s_HodoRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    s_HodoRawHit *q=cont[i];
    if( q->DetectorId()==DetId &&
        q->LayerId()==Layer &&
        q->SegmentId()==Seg ){
      p=q; break;
    }
  }
  if(!p){
    p = new s_HodoRawHit( DetId, Layer, Seg );
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

const s_HodoRHitContainer& s_BeamRawHit::GetsT0RHC() const
{
  return sT0RHC;
}

const s_TrRHitContainer & s_BeamRawHit::GetsBSFTRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxBSFT ) layer=0;
  return sBSFTRHC[layer];
}
