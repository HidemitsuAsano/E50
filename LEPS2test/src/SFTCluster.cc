/*
  SFTCluster.cc

  2016/5 H. Asano
*/

#include "SFTCluster.hh"

#include <cmath>
#include <iostream>
#include <iomanip>

SFTCluster::SFTCluster( ): 
  clustersize_(0),
  clusterlxsize_(0.),
  clusterlzsize_(0),
  clusterID_(-1),
  localx_(-9999.),
  tiltangle_(0.),
  globalz_(-9999.),//,adcsum_(-9999)
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


void SFTCluster::Print( void )
{
   std::cout << "cluster size:" << clustersize_ << std::endl; 
   std::cout << "lx size     :" << clusterlxsize_ << std::endl; 
   std::cout << "lz size     :" << clusterlzsize_ << std::endl; 
   std::cout << "cluster ID  :" << clusterID_ << std::endl;
   std::cout << "local x     :" << localx_    << std::endl;
   std::cout << "tiltangle   :" << tiltangle_ << std::endl;
   std::cout << "global z    :" << globalz_   << std::endl;
   std::cout << "layer       :" << layer_  << std::endl;
   std::cout << "associated? :" << AssociatedLocalTrack_ << std::endl;
   std::cout << "projected X :" << ProjectedPosX_ << std::endl;
   std::cout << "projected Y :" << ProjectedPosY_ << std::endl;
}
