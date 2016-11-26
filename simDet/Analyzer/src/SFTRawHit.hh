/*
  SFTRawHit.hh
  
  2016/2  K.Shirotori
*/

//added comments by H.Asano
//This is the data object class of a single hit for each layer of trackers (SFT and DC ?) 
//actual values is stored in the RawData.cc

#ifndef SFTRawHit_h 
#define SFTRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;

class SFTRawHit
{

public:
  SFTRawHit( int layer, int ch );
  ~SFTRawHit() {};

private:
  int LayerId_, ChId_;
  DoubleVec  PosX_, PosY_, DL_;//probably obsolete, these are only used in simple detector 
  int HitID_;


public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetWire( int ch ) { ChId_=ch; }
  void SetPosX( double posx ) { PosX_.push_back(posx); }
  void SetPosY( double posy ) { PosY_.push_back(posy); }
  void SetDL( double dl ) { DL_.push_back(dl); }
  void SetHitID (int id) { HitID_ = id;};

  int LayerId( void ) const { return LayerId_; }
  int ChId( void ) const { return ChId_; }
  double GetPosX( int nh ) const { return PosX_[nh]; }
  double GetPosY( int nh ) const { return PosY_[nh]; }
  double GetDL( int nh ) const { return DL_[nh]; }
  int GetSize( void ) const { return DL_.size(); }
  int GetHitID( void ) const { return HitID_; } 

};
#endif
