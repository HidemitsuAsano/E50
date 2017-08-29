//This is the data object class of a single hit for each layer of trackers (SFT) 
//actual values are stored in the RawData.cc

// description of member variables:
// LayerId_ : layer number  0-11
// ChId_ : channel number 0-127
// type : X,U or V. X:0,U:1,V:2
// fiberId_ : id of fiber in a layer 
//            in X layer 0-11
//            in U and V layer 0-9
// HitId_ : hit id ? not used now.
// AdcHigh_ : ADC high value
// otrAdcHigh_ : ADC high value is over threshold or not
// AdcLow_ : ADC low value
// otrAdcLow_;  ADC low value is over threshold or not
// TdcLeading_;  TDC leading edge value
// TdcTrailing_; TDC trailing edge value

#ifndef SFTRawHit_h 
#define SFTRawHit_h

#include <cstddef>
//#include <vector>

//comment Aug.22 2017
//deleted vector structure in this class, because in real data, there is no multiple hits in one data.
//typedef std::vector <double> DoubleVec;

class SFTRawHit
{

public:
  SFTRawHit(int ch );
  ~SFTRawHit() {};
private:
  int LayerId_, ChId_;
  int type_;
  int fiberId_;
  //DoubleVec  PosX_, PosY_, DL_;//probably obsolete, these are only used in simple detector 
  int HitID_;
  int AdcHigh_;
  bool otrAdcHigh_;
  int AdcLow_;
  bool otrAdcLow_;
  int TdcLeading_;
  int TdcTrailing_;
 

public:   
  void SetLayer( int layer ) { LayerId_=layer; }
  void SetChId( int ch ) { ChId_=ch; }
  void SetType( int xuv ) { type_ = xuv; }
  void SetFiberId( int id ) { fiberId_ = id; }
  //void SetPosX( double posx ) { PosX_.push_back(posx); }
  //void SetPosY( double posy ) { PosY_.push_back(posy); }
  //void SetDL( double dl ) { DL_.push_back(dl); }
  void SetHitID (int id ) { HitID_ = id; }
  void SetAdcHigh(int val ) { AdcHigh_ = val; }
  void SetotrAdcHigh(bool yes ) { otrAdcHigh_ = yes; }
  void SetAdcLow(int val){AdcLow_ = val;}
  void SetotrAdcLow(bool yes ) { otrAdcLow_ = yes; }
  void SetTdcLeading(int val){TdcLeading_ = val;}
  void SetTdcTrailing(int val){TdcTrailing_ = val;}

  int LayerId( void ) const { return LayerId_; }
  int ChId( void ) const { return ChId_; }
  int FiberId ( void ) const { return fiberId_ ;}
  int GetType ( void ) const { return type_; }
  //double GetPosX( int nh ) const { return PosX_[nh]; }
  //double GetPosY( int nh ) const { return PosY_[nh]; }
  //double GetDL( int nh ) const { return DL_[nh]; }
  //int GetSize( void ) const { return DL_.size(); }
  int GetHitID( void ) const { return HitID_; }
  int GetAdcHigh ( void ) const { return AdcHigh_; }
  bool GetotrAdcHigh ( void ) const { return otrAdcHigh_; }
  int GetAdcLow ( void ) const { return AdcLow_; }
  bool GetotrAdcLow ( void ) const { return otrAdcLow_; }
  int GetTdcLeading ( void ) const { return TdcLeading_; }
  int GetTdcTrailing ( void ) const { return TdcTrailing_; }
  void Print();

};
#endif
