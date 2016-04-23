/*
  s_ScatRawHit.cc

  2016/2
*/

#include "s_ScatRawHit.hh"
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

s_ScatRawHit::~s_ScatRawHit()
{
  clearRegisteredHits();
}

void s_ScatRawHit::clearRegisteredHits( void )
{
  static const std::string funcname="[s_ScatRawHit::clearRegisteredHits]";

  std::for_each(sTOFRHC.begin(), sTOFRHC.end(), DeleteObject());
  sTOFRHC.clear();

  std::for_each(sITOFRHC.begin(), sITOFRHC.end(), DeleteObject());
  sITOFRHC.clear();

  std::for_each(sPADRHC.begin(), sPADRHC.end(), DeleteObject());
  sPADRHC.clear();

  std::for_each(sRICHRHC.begin(), sRICHRHC.end(), DeleteObject());
  sRICHRHC.clear();

  std::for_each(sPID1RHC.begin(), sPID1RHC.end(), DeleteObject());
  sPID1RHC.clear();

  std::for_each(sPID2RHC.begin(), sPID2RHC.end(), DeleteObject());
  sPID2RHC.clear();

  std::for_each(sMFRHC.begin(), sMFRHC.end(), DeleteObject());
  sMFRHC.clear();

  std::for_each(sVDRHC.begin(), sVDRHC.end(), DeleteObject());
  sVDRHC.clear();

  for( int l=0; l<=PlMaxSFT; ++l){
    for_each( sSFTRHC[l].begin(),  sSFTRHC[l].end(), DeleteObject());
    sSFTRHC[l].clear();
  }

  for( int l=0; l<=(PlMaxAFTR+PlMaxAFTL); ++l){
    for_each( sAFTRHC[l].begin(),  sAFTRHC[l].end(), DeleteObject());
    sAFTRHC[l].clear();
  }

  for( int l=0; l<=PlMaxIT1; ++l){
    for_each( sIT1RHC[l].begin(),  sIT1RHC[l].end(), DeleteObject());
    sIT1RHC[l].clear();
  }

  for( int l=0; l<=(PlMaxIT2R+PlMaxIT2L); ++l){
    for_each( sIT2RHC[l].begin(),  sIT2RHC[l].end(), DeleteObject());
    sIT2RHC[l].clear();
  }

  for( int l=0; l<=PlMaxST1; ++l){
    for_each( sST1RHC[l].begin(),  sST1RHC[l].end(), DeleteObject());
    sST1RHC[l].clear();
  }

  for( int l=0; l<=PlMaxST2; ++l){
    for_each( sST2RHC[l].begin(),  sST2RHC[l].end(), DeleteObject());
    sST2RHC[l].clear();
  }
}

bool s_ScatRawHit::SetsTrRHit( int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_ScatRawHit::SetTrRHit]";
  
  //SFT
  if( Layer>=PlMinSFT+PlOffsSFT && Layer<=PlMaxSFT+PlOffsSFT ){
    AddsTrRHit(sSFTRHC[Layer-PlOffsSFT], Layer, Wire, PosX, PosY, DL);
  }
  
  //AFTR
  if( Layer>=PlMinAFTR+PlOffsAFT && Layer<=PlMaxAFTR+PlOffsAFT ){
    AddsTrRHit(sAFTRHC[Layer-PlOffsAFT], Layer, Wire, PosX, PosY, DL);
  }
  
  //AFTL
  if( Layer>=PlMinAFTL+PlOffsAFT+2 && Layer<=PlMaxAFTL+PlOffsAFT+2 ){
    AddsTrRHit(sAFTRHC[Layer-PlOffsAFT], Layer, Wire, PosX, PosY, DL);
  }
  
  //IT1
  if( Layer>=PlMinIT1+PlOffsIT1 && Layer<=PlMaxIT1+PlOffsIT1 ){
    AddsTrRHit(sIT1RHC[Layer-PlOffsIT1], Layer, Wire, PosX, PosY, DL);
  }
  
  //IT2R
  if( Layer>=PlMinIT2R+PlOffsIT2 && Layer<=PlMaxIT2R+PlOffsIT2 ){
    AddsTrRHit(sIT2RHC[Layer-PlOffsIT2], Layer, Wire, PosX, PosY, DL);
  }
  
  //IT2L
  if( Layer>=PlMinIT2L+PlOffsIT2+10 && Layer<=PlMaxIT2R+PlOffsIT2+10 ){
    AddsTrRHit(sIT2RHC[Layer-PlOffsIT2-1], Layer, Wire, PosX, PosY, DL);
  }
  
  //ST1
  if( Layer>=PlMinST1+PlOffsST1 && Layer<=PlMaxST1+PlOffsST1 ){
    AddsTrRHit(sST1RHC[Layer-PlOffsST1], Layer, Wire, PosX, PosY, DL);
  }
  
  //ST2
  if( Layer>=PlMinST2+PlOffsST2 && Layer<=PlMaxST2+PlOffsST2 ){
    AddsTrRHit(sST2RHC[Layer-PlOffsST2], Layer, Wire, PosX, PosY, DL);
  }
  
  return true;
}

bool s_ScatRawHit::SetsHodoRHit( int DetId, int Layer, int Seg,
				 double Time, double Edep, 
				 double Path, double Mom,
				 double PosX, double PosY, 
				 int Pid, double Beta )
{
  static const std::string funcname="[s_ScatRawHit::SetHodoRHit]";

  const TrGeomMan & geomMan=TrGeomMan::GetInstance();

  //TOF
  if( DetId == geomMan.GetDetectorId("TOF") ){
    AddsHodoRHit(sTOFRHC,
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //ITOF
  if( DetId == geomMan.GetDetectorId("ITOFR") ||
      DetId == geomMan.GetDetectorId("ITOFL") ){
    AddsHodoRHit(sITOFRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //ITOFS
  if( DetId ==  geomMan.GetDetectorId("ITOFSR") ||
      DetId ==  geomMan.GetDetectorId("ITOFSL") ){
    AddsHodoRHit(sITOFRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //PAD
  if( DetId == geomMan.GetDetectorId("PADU") ||
      DetId == geomMan.GetDetectorId("PADD") ){
    AddsHodoRHit(sPADRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //RICH
  if( DetId == geomMan.GetDetectorId("RICH") ){
    AddsHodoRHit(sRICHRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //PID1
  if( DetId == geomMan.GetDetectorId("PID1R") ||
      DetId == geomMan.GetDetectorId("PID1L") ){
    AddsHodoRHit(sPID1RHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //PID2
  if( DetId == geomMan.GetDetectorId("PID2R") ||
      DetId == geomMan.GetDetectorId("PID2L") ){
    AddsHodoRHit(sPID2RHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //MF
  if( DetId == geomMan.GetDetectorId("MFU") ||
      DetId == geomMan.GetDetectorId("MFD") ){
    AddsHodoRHit(sMFRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  //VD
  if( DetId == geomMan.GetDetectorId("SpecVp1") ||
      DetId == geomMan.GetDetectorId("SpecVp2") ||
      DetId == geomMan.GetDetectorId("SpecVp3") ){
    AddsHodoRHit(sVDRHC, 
		 DetId, Layer, Seg, Time, Edep, Path, Mom, 
		 PosX, PosY, Pid, Beta);
  }
  
  return true;
}

bool s_ScatRawHit::AddsTrRHit( s_TrRHitContainer &cont,
			       int Layer, int Wire,
			       double PosX, double PosY, double DL)
{
  static const std::string funcname="[s_ScatRawHit::AddTrRHit]";
  
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

bool s_ScatRawHit::AddsHodoRHit( s_HodoRHitContainer &cont,
				 int DetId, int Layer, int Seg,
				 double Time, double Edep, 
				 double Path, double Mom,
				 double PosX, double PosY, 
				 int Pid, double Beta )
{
  static const std::string funcname="[s_ScatRawHit::AddHodoRHit]";
  
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

const s_HodoRHitContainer& s_ScatRawHit::GetsTOFRHC() const
{
  return sTOFRHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsITOFRHC() const
{
  return sITOFRHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsPADRHC() const
{
  return sPADRHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsRICHRHC() const
{
  return sRICHRHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsPID1RHC() const
{
  return sPID1RHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsPID2RHC() const
{
  return sPID2RHC;
}


const s_HodoRHitContainer& s_ScatRawHit::GetsMFRHC() const
{
  return sMFRHC;
}

const s_HodoRHitContainer& s_ScatRawHit::GetsVDRHC() const
{
  return sVDRHC;
}

const s_TrRHitContainer & s_ScatRawHit::GetsSFTRHC( int layer ) const
{
  if( layer<0 || layer>PlMaxSFT ) layer=0;
  return sSFTRHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetsAFTRHC( int layer ) const
{
  if( layer<0 || layer>(PlMaxAFTR+PlMaxAFTL) ) layer=0;
  return sAFTRHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetsIT1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxIT1 ) layer=0;
  return sIT1RHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetsIT2RHC( int layer ) const
{
  if( layer<0 || layer>(PlMaxIT2R+PlMaxIT2L) ) layer=0;
  return sIT2RHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetsST1RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxST1 ) layer=0;
  return sST1RHC[layer];
}

const s_TrRHitContainer & s_ScatRawHit::GetsST2RHC( int layer ) const
{
  if( layer<0 || layer>PlMaxST2 ) layer=0;
  return sST2RHC[layer];
}
