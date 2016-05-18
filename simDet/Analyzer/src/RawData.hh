/*
  RawData.hh

  2016/4  K.Shirotori
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

typedef std::vector<PrimInfo*>   PrimInfoContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<TrRawHit*>   TrRHitContainer;

class RawData
{

private:
  PrimInfoContainer PrimHC;

  HodoRHitContainer T0RHC;

  TrRHitContainer SFTRawHitContainer[PlMaxSFT+1];//vector of data object class for a single track in each SFT layer, PlMaxSFT: max layer number of SFT (=12)

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

public:
  const PrimInfoContainer& GetPrimHC() const;

  const HodoRHitContainer& GetT0RHC() const;

  const TrRHitContainer & GetSFTRawHitContainer( int layer ) const;

};

#endif

