/*
  TrHitCluster.cc

  2012/5 K.Shirotori
*/

#include "TrHitCluster.hh"

TrHitCluster::TrHitCluster( TrLTrackHit *hitA, TrLTrackHit *hitB )
  : hitA_(hitA), hitB_(hitB), clustersize_(0)
{
  if(hitA_) ++clustersize_;
  if(hitB_) ++clustersize_;
}

TrHitCluster::~TrHitCluster()
{
}
