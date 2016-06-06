/*
  RawData.hh

  2016/2  K.Shirotori
*/

#ifndef RawData_h
#define RawData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DetectorID.hh"
#include <vector>

class PrimInfo;
class HodoRawHit;
class TrRawHit;

class s_BeamRawHit;
class s_ScatRawHit;

typedef std::vector<PrimInfo*>   PrimInfoContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<TrRawHit*>   TrRHitContainer;

typedef std::vector<s_BeamRawHit*> s_BeamRHitContainer;
typedef std::vector<s_ScatRawHit*> s_ScatRHitContainer;

class RawData
{

private:
  PrimInfoContainer PrimHC;

  HodoRHitContainer T0RHC;
  HodoRHitContainer TOFRHC;
  HodoRHitContainer ITOFRHC;
  HodoRHitContainer PADRHC;
  HodoRHitContainer RICHRHC;
  HodoRHitContainer PID1RHC;
  HodoRHitContainer PID2RHC;
  HodoRHitContainer MFRHC;
  HodoRHitContainer VDRHC;

  TrRHitContainer BFTRHC[PlMaxBFT+1];
  TrRHitContainer SFTRHC[PlMaxSFT+1];
  TrRHitContainer AFTRHC[PlMaxAFTR+PlMaxAFTL+1];
  TrRHitContainer IT1RHC[PlMaxIT1+1];
  TrRHitContainer IT2RHC[PlMaxIT2R+PlMaxIT2L+1];
  TrRHitContainer ST1RHC[PlMaxST1+1];
  TrRHitContainer ST2RHC[PlMaxST2+1];

  s_BeamRHitContainer s_BeamRHC;
  s_ScatRHitContainer s_ScatRHC;

public:
  RawData();
  ~RawData();

  void clearAll();
  bool DecodeRawHits( std::ifstream & );

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  //PrimInfo
  bool AddPrimInfo( PrimInfoContainer& cont,
		    double VertX, double VertY, double VertZ,
		    double BeamMom, double BeamU, double BeamV,
		    double AnaBeamMom,
		    double Mass1, double Mom1,
		    double Theta1, double Phi1, double ThetaCM1, double PhiCM1,
		    double Mass2, double Mom2,
		    double Theta2, double Phi2, double ThetaCM2, double PhiCM2 );

  //Full tracking
  bool AddTrRHit( TrRHitContainer &cont, 
		  int Layer, int Wire, 
		  double PosX, double PosY, double DL );
  
  bool AddHodoRHit( HodoRHitContainer& cont,
		    int DetId, int Layer, int Seg,
		    double Time, double Edep, 
		    double Path, double Mom,
		    double PosX, double PosY, 
		    int Pid, double Beta );

  //Simple tracking
  bool AddsBGeneRHit( s_BeamRHitContainer& cont, 
		      int TrackID, int Type,
		      double Mom, 
		      double MomX, double MomY, double MomZ,  
		      int Pid, 
		      double PosX, double PosY );

  bool AddsBTrRHit( s_BeamRHitContainer& cont, 
		    int TrackID, int Type,
		    int Layer, int Wire, 
		    double PosX, double PosY, double DL );
  
  bool AddsBHodoRHit( s_BeamRHitContainer& cont,
		      int TrackID, int Type, 
		      int DetId, int Layer, int Seg,
		      double Time, double Edep, 
		      double Path, double Mom,
		      double PosX, double PosY, 
		      int Pid, double Beta );

  bool AddsSGeneRHit( s_ScatRHitContainer& cont, 
		      int TrackID, int Type, 
		      double Mom, 
		      double MomX, double MomY, double MomZ,  
		      int Pid, 
		      double PosX, double PosY );
  
  bool AddsSTrRHit( s_ScatRHitContainer& cont, 
		    int TrackID, int Type, 
		    int Layer, int Wire, 
		    double PosX, double PosY, double DL );

  bool AddsSHodoRHit( s_ScatRHitContainer& cont,
		      int TrackID, int Type,
		      int DetId, int Layer, int Seg,
		      double Time, double Edep, 
		      double Path, double Mom,
		      double PosX, double PosY, 
		      int Pid, double Beta );

public:
  const PrimInfoContainer& GetPrimHC() const;

  const HodoRHitContainer& GetT0RHC() const;
  const HodoRHitContainer& GetTOFRHC() const;
  const HodoRHitContainer& GetITOFRHC() const;
  const HodoRHitContainer& GetPADRHC() const;
  const HodoRHitContainer& GetRICHRHC() const;
  const HodoRHitContainer& GetPID1RHC() const;
  const HodoRHitContainer& GetPID2RHC() const;
  const HodoRHitContainer& GetMFRHC() const;
  const HodoRHitContainer& GetVDRHC() const;

  const TrRHitContainer & GetBFTRHC( int layer ) const;
  const TrRHitContainer & GetSFTRHC( int layer ) const;
  const TrRHitContainer & GetAFTRHC( int layer ) const;
  const TrRHitContainer & GetIT1RHC( int layer ) const;
  const TrRHitContainer & GetIT2RHC( int layer ) const;
  const TrRHitContainer & GetST1RHC( int layer ) const;
  const TrRHitContainer & GetST2RHC( int layer ) const;

  const s_BeamRHitContainer & GetsBeamRHC() const;
  const s_ScatRHitContainer & GetsScatRHC() const;

};

#endif

