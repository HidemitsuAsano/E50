/*
  HodoAnalyzer.cc

  2016/2  K.Shirotori
*/

#include "HodoAnalyzer.hh"
#include "RawData.hh"
#include "HodoHit.hh"
#include "HodoCluster.hh"

#include "TemplateLib.hh"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

const double MaxTimeDifT0  = 0.1;
const double MaxTimeDifTOF = 0.1;
const double MaxTimeDifITOF= 0.1;
const double MaxTimeDifPAD = 0.1;

#define Cluster 1

HodoAnalyzer::HodoAnalyzer()
{
}

HodoAnalyzer::~HodoAnalyzer()
{
  clearT0Hits();
}

void HodoAnalyzer::clearT0Hits()
{
  for_each(T0Cont.begin(),T0Cont.end(),DeleteObject());
  for_each(T0ClCont.begin(),T0ClCont.end(),DeleteObject());
}

bool HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeT0Hits( rawData );

  return true;
}

bool HodoAnalyzer::DecodeT0Hits( RawData *rawData )
{
  clearT0Hits();
  
  const HodoRHitContainer &cont=rawData->GetT0RHC();
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit=cont[i];
    if( !hit ) continue;
    HodoHit *hp=new HodoHit( hit );
    if( !hp ) continue;
    if( hp->calculate() )
      T0Cont.push_back(hp);
    else
      delete hp;
  }
  
#if Cluster
  MakeUpClusters( T0Cont, T0ClCont, MaxTimeDifT0 );
#endif
  
  return true;
}

int HodoAnalyzer::
MakeUpClusters( const HodoHitContainer & HitCont,
		HodoClusterContainer & ClusterCont, double maxTimeDif )
{
  static const std::string funcname = "HodoAnalyzer::MakeUpClusters";
  
  if( !ClusterCont.empty() )
    for_each( ClusterCont.begin(), ClusterCont.end(), DeleteObject() );
  
  int nh=HitCont.size();
  
  std::vector <int> flag(nh,0);
  
  for(int i=0; i<nh; ++i ){
    if( flag[i] ) continue;
    HodoHit *hitA=HitCont[i];
    int segA=hitA->SegmentId();
    double cmtA=hitA->GetTime(0);
    HodoHit *hitB=0;
    int iB=-1;
    double cmtB;
    int segB;
    for( int j=i+1; j<nh; ++j ){
      HodoHit *hit=HitCont[j];
      int seg=hit->SegmentId();
      double cmt=hit->GetTime(0);
      if( abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif ){
	hitB=hit; ++flag[j]; iB=j; segB=seg; cmtB=cmt; break;
      }
    }
    if(hitB){
      HodoHit *hitC=0;
      for( int j=i+1; j<nh; ++j ){
	if( j==iB ) continue;
	HodoHit *hit=HitCont[j];
	int seg=hit->SegmentId();
	double cmt=hit->GetTime(0);
	if( (abs(seg-segA)==1 && fabs(cmt-cmtA)<maxTimeDif) ||
	    (abs(seg-segB)==1 && fabs(cmt-cmtB)<maxTimeDif) ){
	  hitC=hit; ++flag[j]; break;
	}
      }
      if(hitC){
	HodoCluster *cluster=new HodoCluster(hitA,hitB,hitC);
	if( cluster ) ClusterCont.push_back(cluster);
      }
      else{
	HodoCluster *cluster=new HodoCluster(hitA,hitB);
	if( cluster ) ClusterCont.push_back(cluster);
      }
    }
    else{
      HodoCluster *cluster=new HodoCluster(hitA);
      if( cluster ) ClusterCont.push_back(cluster);
    }
  }

  return ClusterCont.size(); 
}				 

bool HodoAnalyzer::ReCalcT0Hits( bool applyRecursively )
{
  int n=T0Cont.size();
  for( int i=0; i<n; ++i ){
    HodoHit *hit=T0Cont[i];
    if(hit) hit->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcT0Clusters( bool applyRecursively )
{
  int n=T0ClCont.size();
  for( int i=0; i<n; ++i ){
    HodoCluster *cl=T0ClCont[i];
    if(cl) cl->ReCalc(applyRecursively);
  }
  return true;
}

bool HodoAnalyzer::ReCalcAll( void )
{
  ReCalcT0Hits();

  ReCalcT0Clusters();  

  return true;
}

