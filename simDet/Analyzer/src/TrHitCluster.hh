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
  float clusterlxsize_; // cluster size in local x axis
  int clusterlzsize_; // 1 or 2 (max number = number of sublayer)
  unsigned int clusterID_;
  float localx_;//mean (maybe weighted mean) of TrHits
//  float adcsum_;//sum of adc of each MPPC
                 //TODO : how to deal with TDC info ?
  float tiltangle_; // maybe don't need to implement here, since tiltangle and z position can be obtained from the confmanager
  float localz_;      //

public:
  unsigned int GetClusterID(void) const { return clusterID_; }
  void SetClusterID(unsigned int id)  { clusterID_ = id; }
  int GetClusterSize( void ) const { return clustersize_; }
  void SetClusterSize(int size) { clustersize_ = size; }
  float GetClusterLxSize( void ) const { return clusterlxsize_; }
  void SetClusterLxSize(float lxsize) { clusterlxsize_ = lxsize; }
  int GetClusterLzSize( void ) const { return clusterlzsize_; }
  void SetClusterLzSize(int lzsize) { clusterlzsize_ = lzsize; }
  float GetLocalX( void ) const { return localx_; }
  void SetLocalX(float lx )  { localx_=lx; }
  float GetLocalZ( void ) const { return localz_; }
  void SetLocalZ(float lz )  { localz_=lz; }
  void SetTiltAngle (float angle) { tiltangle_ = angle; }
  float GetTiltAngle (void) { return tiltangle_ ; }
//  float GetAdcSum( void ) const { return adcsum_; }
//  void SetAdcSum(float adc )  { adcsum_=adc; }
  /*
  TrLTrackHit *GetHit( int i ) const
  {
    if(i==0)      return hitA_;
    else if(i==1) return hitB_;
    else          return 0;
  }*/

};

#endif
