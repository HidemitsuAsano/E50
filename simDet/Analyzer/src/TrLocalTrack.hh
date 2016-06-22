/*
  TrLocalTrack.hh

  2012/5  K.Shirotori
*/

#ifndef TrLocalTrack_h
#define TrLocalTrack_h 1

#include <vector>
#include <functional>
#include "ThreeVector.hh"

class TrLTrackHit;
class TrAnalyzer;

class TrLocalTrack
{
public:
  explicit TrLocalTrack();
  ~TrLocalTrack();
private:
  TrLocalTrack( const TrLocalTrack & );
  TrLocalTrack & operator = ( const TrLocalTrack & );

private:
  std::vector <TrLTrackHit *> hitArray;
  /*
  double Av_;//not used
  double Ax_;//not used
  double Au_;//not used
  double Chiv_;//not used
  double Chix_;//not used
  double Chiu_;//not used
  */
public:
  void AddHit( TrLTrackHit *hitp ) { hitArray.push_back( hitp ); }
  bool DoFit( void );//Liner func
  bool DoFit2( void );//Quadrotic func
  std::size_t GetNHit( void ) const { return hitArray.size(); }
  TrLTrackHit * GetHit( std::size_t nth ) const;
  TrLTrackHit * GetHitOfLayerNumber( int lnum ) const;

  // optional extension by Yone
  bool DoFitVXU( void );//Liner func
  bool DoFitVXU2( void );//Quadrotic func
  //void SetAv( double Av) { Av_=Av; }
  //void SetAx( double Ax) { Ax_=Ax; }
  //void SetAu( double Au) { Au_=Au; }
  //void SetChiv( double Chiv) { Chiv_=Chiv; }
  //void SetChix( double Chix) { Chix_=Chix; }
  //void SetChiu( double Chiu) { Chiu_=Chiu; }

  double GetX0( void ) const { return x0_; }
  double GetY0( void ) const { return y0_; }
  double GetU0( void ) const { return u0_; }
  double GetV0( void ) const { return v0_; }
  double GetU1( void ) const { return u1_; }
  double GetV1( void ) const { return v1_; }

  double GetVXU_A( void ) const { return a_; }
  double GetVXU_B( void ) const { return b_; }
  double GetVXU_C( void ) const { return c_; }
  //double GetAv( void ) const { return Av_; }
  //double GetAx( void ) const { return Ax_; }
  //double GetAu( void ) const { return Au_; }
  
  /*
  double GetDifVXU( void ) 
    const { return (Av_/cos(acos(-1.)/180.*(-15.0))-Ax_)*(Av_/cos(acos(-1.)/180.*(-15.0))-Ax_)
      +(Ax_-Au_/cos(acos(-1.)/180.*(15.0)))*(Ax_-Au_/cos(acos(-1.)/180.*(15.0)))
      +(Au_/cos(acos(-1.)/180.*(15.0))
	-Av_/cos(acos(-1.)/180.*(-15.0)))*(Au_/cos(acos(-1.)/180.*(15.0))
					   -Av_/cos(acos(-1.)/180.*(-15.0)) ) ; }

  double GetDifVXU2( void ) 
    const { return (Av_*cos(acos(-1.)/180.*(-30.0))-Ax_)*(Av_*cos(acos(-1.)/180.*(-30.0))-Ax_)
      +(Ax_-Au_*cos(acos(-1.)/180.*(30.0)))*(Ax_-Au_*cos(acos(-1.)/180.*(30.0)))
      +(Au_*cos(acos(-1.)/180.*(30.0))
	-Av_*cos(acos(-1.)/180.*(-30.0)))*(Au_*cos(acos(-1.)/180.*(30.0))
					   -Av_*cos(acos(-1.)/180.*(-30.0)))
      +(Au_*sin(acos(-1.)/180.*(30.0))
	-Av_*sin(acos(-1.)/180.*(30.0)))*(Au_*sin(acos(-1.)/180.*(30.0))
					  -Av_*sin(acos(-1.)/180.*(30.0))) ; }
  */
  double GetChiSquare( void ) const { return chisqr_; }
  //double GetChiV( void ) const { return Chiv_; }
  //double GetChiX( void ) const { return Chix_; }
  //double GetChiU( void ) const { return Chiu_; }
  double GetX( double z ) const { return x0_+u0_*z+u1_*z*z; } 
  double GetY( double z ) const { return y0_+v0_*z+v1_*z*z; } 
  bool GetStatus( void ) const { return status_; } 
  bool GoodForTracking( void ) const { return gftstatus_; }
  bool GoodForTracking( bool status )
  { bool ret=gftstatus_; gftstatus_=status; return ret; } 
  bool ReCalc( bool ApplyRecursively=false );  

private:
  bool status_;
  double x0_, y0_, u0_, v0_, u1_, v1_; // x0_, y0_ .. : fitting results of beam position at (0,0,0) in SFT coordinate system.
  double a_,b_,c_;
  double chisqr_;
  bool gftstatus_;

  double *Coefficients_x[6];
  double *Coefficients_y[6];
  inline bool allocateBufferArea( void );

};

struct TrLTrackComp 
  : public std::binary_function <TrLocalTrack *, TrLocalTrack *, bool>
{
  bool operator()( const TrLocalTrack * const p1, 
		   const TrLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    if( (n1>n2+1) ){
      return true;
    }
    else if( (n2>n1+1)  ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};

struct TrLTrackComp1 
  : public std::binary_function <TrLocalTrack *, TrLocalTrack *, bool>
{
  bool operator()( const TrLocalTrack * const p1, 
		   const TrLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    if(n1>n2) return true;
    else if(n2>n1) return false;
    else
      return (p1->GetChiSquare())<=(p2->GetChiSquare());
  }

};

struct TrLTrackComp2 
  : public std::binary_function <TrLocalTrack *, TrLocalTrack *, bool>
{
  bool operator()( const TrLocalTrack * const p1, 
		   const TrLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    if(n1<n2) return true;
    else if(n2<n1) return false;
    else
      return (p1->GetChiSquare())<=(p2->GetChiSquare());
  }

};

struct TrLTrackComp3 
  : public std::binary_function <TrLocalTrack *, TrLocalTrack *, bool>
{
  bool operator()( const TrLocalTrack * const p1, 
		   const TrLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    double a1=fabs(1.-chi1),a2=fabs(1.-chi2);
    if(a1<a2) return true;
    else if(a2<a1) return false;
    else
      return (n1<=n2);
  }

};

struct TrLTrackComp4 
  : public std::binary_function <TrLocalTrack *, TrLocalTrack *, bool>
{
  bool operator()( const TrLocalTrack * const p1, 
		   const TrLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();
    //if( (n1>n2+1) ){
    //    if( (n1>n2+1) && (fabs(chi1-chi2)<5.) ){
    if( (n1>n2+1) && (fabs(chi1-chi2)<2.) ){
      return true;
    }
    //else if( (n2>n1+1)  ){
    //    else if( (n2>n1+1) && (fabs(chi1-chi2)<5.) ){
    else if( (n2>n1+1) && (fabs(chi1-chi2)<2.) ){
      return false;
    }
    else{
      return (chi1<=chi2);
    }
  }
};



#endif
