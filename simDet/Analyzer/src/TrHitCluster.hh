/*
  TrHitCluster.hh

  2012/5 K.Shirotori
*/

#ifndef TrHitCluster_h
#define TrHitCluster_h 1

//class TrLTrackHit;

class TrHitCluster
{
public:
  //TrHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB=0 );
  TrHitCluster( );
  ~TrHitCluster();

private:
  //TrLTrackHit *hitA_, *hitB_;
  int clustersize_; // number of hit fiber (MPPC)
  int clusterlxsize_; // cluster size in local x axis
  int clusterlzsize_; // 1 or 2 (max number = number of sublayer)
  unsigned int clusterID_;
  double localx_;//mean (maybe weighted mean) of TrHits
  double adcsum_;//sum of adc of each MPPC
                 //TODO : how to deal with TDC info ?
  double tiltangle_; // maybe don't need to implement here, since tiltangle and z position can be obtained from the confmanager
  double zpos_;      //

public:
  unsigned int GetClusterID(void) const { return clusterID_; }
  void SetClusterID(unsigned int id)  { clusterID_ = id; }
  int GetClusterSize( void ) const { return clustersize_; }
  void SetClusterSize(int size) { clustersize_ = size; }
  int GetClusterLxSize( void ) const { return clusterlxsize_; }
  void SetClusterLxSize(int lxsize) { clusterlxsize_ = lxsize; }
  int GetClusterLzSize( void ) const { return clusterlzsize_; }
  void SetClusterLzSize(int lzsize) { clusterlzsize_ = lzsize; }
  double GetLocalX( void ) const { return localx_; }
  void SetLocalX(double lx )  { localx_=lx; }
  double GetAdcSum( void ) const { return adcsum_; }
  void SetAdcSum(double adc )  { adcsum_=adc; }
  
  /*
  TrLTrackHit *GetHit( int i ) const
  {
    if(i==0)      return hitA_;
    else if(i==1) return hitB_;
    else          return 0;
  }*/

};

#endif
