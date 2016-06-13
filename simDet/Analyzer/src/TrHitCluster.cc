/*
  TrHitCluster.cc

  2012/5 K.Shirotori
*/

#include "TrHitCluster.hh"

TrHitCluster::TrHitCluster( ): 
  clustersize_(0),
  clusterlxsize_(0),
  clusterlzsize_(0),
  clusterID_(0),
  localx_(-9999),
  tiltangle_(0),
  localz_(-9999)//,adcsum_(-9999)
{
}

TrHitCluster::~TrHitCluster()
{
}
