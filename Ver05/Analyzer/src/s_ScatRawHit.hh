/*
  s_ScatRawHit.hh

  2016/2
*/

#ifndef s_ScatRawHit_h 
#define s_ScatRawHit_h 1

#include "DetectorID.hh"
#include <vector>

class s_TrRawHit;
class s_HodoRawHit;

typedef std::vector<s_TrRawHit*>   s_TrRHitContainer;
typedef std::vector<s_HodoRawHit*> s_HodoRHitContainer;

class s_ScatRawHit
{
private:
  int TrackId_, TrackType_;

  s_HodoRHitContainer sTOFRHC;
  s_HodoRHitContainer sITOFRHC;
  s_HodoRHitContainer sPADRHC;
  s_HodoRHitContainer sRICHRHC;
  s_HodoRHitContainer sPID1RHC;
  s_HodoRHitContainer sPID2RHC;
  s_HodoRHitContainer sMFRHC;
  s_HodoRHitContainer sVDRHC;

  s_TrRHitContainer sSSFTRHC[PlMaxSSFT+1];
  s_TrRHitContainer sAFTRHC[PlMaxAFTR+PlMaxAFTL+1];
  s_TrRHitContainer sIT1RHC[PlMaxIT1+1];
  s_TrRHitContainer sIT2RHC[PlMaxIT2R+PlMaxIT2L+1];
  s_TrRHitContainer sST1RHC[PlMaxST1+1];
  s_TrRHitContainer sST2RHC[PlMaxST2+1];

  //Generated
  double p_, px_, py_, pz_;
  int pid_;
  double x_, y_;  

public:
  s_ScatRawHit( int trackid, int tracktype )
    : TrackId_(trackid), TrackType_(tracktype),
      p_(0.0), px_(0.0), py_(0.0), pz_(0.0),
      pid_(0),
      x_(0.0), y_(0.0)
  {}
  ~s_ScatRawHit();

public:
  int TrackId( void ) const { return TrackId_; }
  int TrackType( void ) const { return TrackType_; }

  //Generated
  void SetMom( double  p ) { p_=p; }
  void SetMomX( double px ) { px_=px; }
  void SetMomY( double py ) { py_=py; }
  void SetMomZ( double pz ) { pz_=pz; }
  void SetPid( int pid ) { pid_=pid; }
  void SetVertX( double x ) { x_=x; }
  void SetVertY( double y ) { y_=y; }

  double GetMom( void ) const { return p_; };
  double GetMomX( void ) const { return px_; };
  double GetMomY( void ) const { return py_; };
  double GetMomZ( void ) const { return pz_; };
  int GetPid( void ) const { return pid_; };
  double GetVertX( void ) const { return x_; };
  double GetVertY( void ) const { return y_; };

  //Beam information
  const s_HodoRHitContainer& GetsTOFRHC() const;
  const s_HodoRHitContainer& GetsITOFRHC() const;
  const s_HodoRHitContainer& GetsPADRHC() const;
  const s_HodoRHitContainer& GetsRICHRHC() const;
  const s_HodoRHitContainer& GetsPID1RHC() const;
  const s_HodoRHitContainer& GetsPID2RHC() const;
  const s_HodoRHitContainer& GetsMFRHC() const;
  const s_HodoRHitContainer& GetsVDRHC() const;

  const s_TrRHitContainer & GetsSSFTRHC( int layer ) const;
  const s_TrRHitContainer & GetsAFTRHC( int layer ) const;
  const s_TrRHitContainer & GetsIT1RHC( int layer ) const;
  const s_TrRHitContainer & GetsIT2RHC( int layer ) const;
  const s_TrRHitContainer & GetsST1RHC( int layer ) const;
  const s_TrRHitContainer & GetsST2RHC( int layer ) const;

  bool SetsTrRHit( int Layer, int Wire, 
		   double PosX, double PosY, double DL );
  
  bool SetsHodoRHit( int DetId, int Layer, int Seg,
		     double Time, double Edep, 
		     double Path, double Mom,
		     double PosX, double PosY, 
		     int Pid, double Beta );

  bool AddsTrRHit( s_TrRHitContainer& cont,
		   int Layer, int Wire, 
		   double PosX, double PosY, double DL );
  
  bool AddsHodoRHit( s_HodoRHitContainer& cont,
		     int DetId, int Layer, int Seg,
		     double Time, double Edep, 
		     double Path, double Mom,
		     double PosX, double PosY, 
		     int Pid, double Beta );
  
  void clearRegisteredHits( void );

};
  
#endif
