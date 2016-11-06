/*
  HodoRawHit.hh
  
  2016/2  K.Shirotori
*/

#ifndef HodoRawHit_h 
#define HodoRawHit_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <int> IntVec;

class HodoRawHit
{
private:
  int    DetId_, LayerId_, SegId_;
  DoubleVec Time_, Edep_;
  DoubleVec Path_, Mom_;
  DoubleVec PosX_, PosY_;
  IntVec Pid_;
  DoubleVec Beta_;

public:
  HodoRawHit( int detid, int layerid, int segid )
    : DetId_(detid), LayerId_(layerid), SegId_(segid),
      Time_(0.0), Edep_(0.0),
      Path_(0.0), Mom_(0.0),
      PosX_(0.0), PosY_(0.0), 
      Pid_(0), Beta_(0.0)
  {};
  ~HodoRawHit() {};

public:
  void SetTime( double time ) { Time_.push_back(time); }
  void SetEdep( double edep ) { Edep_.push_back(edep); }
  void SetPath( double path ) { Path_.push_back(path); }
  void SetMom( double mom ) { Mom_.push_back(mom); }
  void SetPosX( double posx ) { PosX_.push_back(posx); }
  void SetPosY( double posy ) { PosY_.push_back(posy); }
  void SetPid( int pid ) { Pid_.push_back(pid); }
  void SetBeta( double beta ) { Beta_.push_back(beta); }

  int DetectorId( void ) const { return DetId_; };
  int LayerId( void ) const { return LayerId_; };
  int SegmentId( void ) const { return SegId_; };

  int    GetSize( void ) const { return Time_.size(); };
  double GetTime( int nh ) const { return Time_[nh]; };
  double GetEdep( int nh ) const { return Edep_[nh]; };
  double GetPath( int nh ) const { return Path_[nh]; };
  double GetMom( int nh ) const { return Mom_[nh]; };
  double GetPosX( int nh ) const { return PosX_[nh]; };
  double GetPosY( int nh ) const { return PosY_[nh]; };
  int    GetPid( int nh ) const { return Pid_[nh]; };
  double GetBeta( int nh ) const { return Beta_[nh]; };
};
#endif
