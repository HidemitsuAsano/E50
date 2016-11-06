/*
  RawHit.hh
  
  2016/2  K.Shirotori
*/

//added comments by H.Asano
//This is the data object class of a single hit for each layer of trackers (SFT and DC ?) 
//actual values is stored in the RawData.cc

#ifndef RawHit_h 
#define RawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <bool> BoolVec;

class TrLTrackHit;

class RawHit
{
private:
  int LayerId_, WireId_;
  DoubleVec  PosX_, PosY_, DL_;

public:
  RawHit( int layer, int wire )
    : LayerId_(layer), WireId_(wire), 
      PosX_(0.0), PosY_(0.0), DL_(0.0)
  {};
  ~RawHit() {};

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetWire( int wire ) { WireId_=wire; }
  void SetPosX( double posx ) { PosX_.push_back(posx); }
  void SetPosY( double posy ) { PosY_.push_back(posy); }
  void SetDL( double dl ) { DL_.push_back(dl); }

  int LayerId( void ) const { return LayerId_; }
  int WireId( void ) const { return WireId_; }
  double GetPosX( int nh ) const { return PosX_[nh]; }
  double GetPosY( int nh ) const { return PosY_[nh]; }
  double GetDL( int nh ) const { return DL_[nh]; }
  int GetSize( void ) const { return DL_.size(); }

};
#endif
