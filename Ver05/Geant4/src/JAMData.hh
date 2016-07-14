/*
  JAMData.hh
  
  2012/10  K.Shirotori
*/

#ifndef JAMData_h 
#define JAMData_h

#include <cstddef>
#include <vector>

typedef std::vector <double> DoubleVec;
typedef std::vector <int> IntVec;

class JAMData
{
private:
  int EvNum_;
  IntVec Pid_;
  DoubleVec MomX_, MomY_, MomZ_;

public:
  JAMData( int evnum )
    : EvNum_(evnum), Pid_(0), 
      MomX_(0), MomY_(0), MomZ_(0) 
  {};
  ~JAMData() {};

public:
  void SetPid( int pid ) { Pid_.push_back(pid); }
  void SetMomX( double momx ) { MomX_.push_back(momx); }
  void SetMomY( double momy ) { MomY_.push_back(momy); }
  void SetMomZ( double momz ) { MomZ_.push_back(momz); }

  int EventNo( void ) const { return EvNum_; };

  int    GetPid( int nh ) const { return Pid_[nh]; };
  int    GetPidSize( void ) const { return Pid_.size(); };
  double GetMomX( int nh ) const { return MomX_[nh]; };
  double GetMomY( int nh ) const { return MomY_[nh]; };
  double GetMomZ( int nh ) const { return MomZ_[nh]; };
};
#endif
