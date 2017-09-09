/*
  SFTCluster.hh

  2015/6 H. Asano
*/

#ifndef SFTCluster_h
#define SFTCluster_h 1


class SFTCluster
{
public:
  SFTCluster( );
  ~SFTCluster();

private:
  int clustersize_; // number of hit fibers (or bundles of fibers)
  float clusterlxsize_; // cluster size in local x axis
  int clusterlzsize_; // 1 or 2 (max number = number of sublayer)
  unsigned int clusterID_; // assigned layer by layer , start from 0
  float localx_;//center position of the cluster. The coordinate system is defined in each layer
  // float adcsum_;//sum of adc of each MPPC
                 //TODO : how to deal with TDC info ?
  float tiltangle_; // maybe don't need to implement here, since tiltangle and z position can be obtained from the confmanager (unit : degree)
  float globalz_;    //
  int layer_;
  int AssociatedLocalTrack_;//pointer to the associated local track (TrLocalTrack);
                            //if -1 , there is no associated track
  double ProjectedPosX_;//projected position from tracking x in global coordinate
  double ProjectedPosY_;//projected position from tracking y in global coordinate
  
  
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
  float GetGlobalZ( void ) const { return globalz_; }
  void SetGlobalZ(float lz )  { globalz_=lz; }
  float GetTiltAngle (void) const { return tiltangle_ ; }
  void SetTiltAngle (float angle) { tiltangle_ = angle; }
  int GetLayer(void) const { return layer_;}
  void SetLayer(int lay) { layer_ = lay;}
  int GetAssociatedLocalTrack (void ) const { return AssociatedLocalTrack_ ;}
  void SetAssociatedLocalTrack (int trk = 0) { AssociatedLocalTrack_ = trk;}
  double GetProjectedPosX( void ) const {return ProjectedPosX_;}
  double GetProjectedPosY( void ) const {return ProjectedPosY_;}
  void SetProjectedPosition(double xpos,double ypos)  {ProjectedPosX_ = xpos; ProjectedPosY_ = ypos;}
  double GetProjectedLocalX( void ) const;
  double GetResidual( void ) const { return localx_ - GetProjectedLocalX() ; }
  //  float GetAdcSum( void ) const { return adcsum_; }
  //  void SetAdcSum(float adc )  { adcsum_=adc; }
  void Print( void );

};

#endif
