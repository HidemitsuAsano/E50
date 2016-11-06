/*
  SFTCluster.cc

  2016/5 H. Asano
*/

#include "SFTCluster.hh"

#include <cmath>


SFTCluster::SFTCluster( ): 
  clustersize_(0),
  clusterlxsize_(0.),
  clusterlzsize_(0),
  clusterID_(-1),
  localx_(-9999.),
  tiltangle_(0.),
  localz_(-9999.),//,adcsum_(-9999)
  layer_(-1),
  AssociatedLocalTrack_(-1),
  ProjectedPosX_(-9999.),
  ProjectedPosY_(-9999.)
{
}

SFTCluster::~SFTCluster()
{
}


double SFTCluster::GetProjectedLocalX( void ) const
{
  const double Deg2Rad = acos(-1)/180.;
  double prolocalx = 
  ProjectedPosX_*cos(tiltangle_*Deg2Rad) + ProjectedPosX_*sin(tiltangle_*Deg2Rad);//need degree-> radian

  return prolocalx ;
}
