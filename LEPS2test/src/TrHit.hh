#ifndef TrHit_h 
#define TrHit_h
//TrHit.hh
//This class stores raw hits which have ADC and TDC values above threshold.


//#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <functional>

//comment Aug.22 2017
//deleted vector structure in this class, because in real data, there is no multiple hits in one data.
//typedef std::vector <bool> BoolVec;
//typedef std::vector <int> IntVec;
//typedef std::vector <double> DoubleVec;

//class TrLTrackHit;

class TrHit
{
public:
  TrHit();
  TrHit( int layer, int fiber ); 

  ~TrHit();

private:
  TrHit( const TrHit & );
  TrHit & operator = ( const TrHit & );

private:
  int layer_;
  int fiber_;//fiber id 0-11
  int hitcounter_;
  int adc_;
  int tot_;

  double localx_; //local x position
  double angle_;  //tilt angle

public:
  //bool CalcObservables( void );

  void SetLayer( int layer ) { layer_=layer; }
  int GetLayer( void ) const { return layer_; }
  void SetFiber( int fiber ) { fiber_=fiber; }  
  int GetFiber( void )  const { return fiber_; }
  
  void SetAdc( int val ) { adc_= val;}
  int GetAdc(void) const { return adc_;}

  void SetTot( int val ) { tot_= val;}
  int GetTot( void ) const { return tot_;}
  
  void SetTiltAngle( double angleDegree ) { angle_=angleDegree; }
  double GetTiltAngle( void ) const { return angle_; }
  
  void SetFiberPosition( double pos ) { localx_ = pos; }
  double GetFiberPosition( void ) const { return localx_; }
  
  void Print( void );


  
  static bool compareTrHitPredicate(TrHit *lhs, TrHit *rhs) { return (lhs->fiber_ < rhs->fiber_); }


private:
  void clearRegisteredHits( void );

};


#endif
