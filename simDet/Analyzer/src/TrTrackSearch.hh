/*
  TrTrackSearch.hh

  2012/5  K.Shirotori
*/

#ifndef TrTrackSearch_h
#define TrTrackSearch_h 1

#include "TrAnalyzer.hh"

#include <vector>

class SFTCluster;
class TrLocalTrack;
class TrLTrackHit;

//Full combination by Linear fitting
//int LocalTrackSearch( const TrHitContainer * HC,
int LocalTrackSearch(const  SFTClusterContainer  *HC,
		      std::vector <TrLocalTrack *> &TrackCont,
		      int NumOfLayers, unsigned int MinNumOfHits);

std::vector< std::vector<int> > makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_below( int ndim_org, int maximumHit, int ndim, const int *index1 ); 


TrLocalTrack *MakeTrack(const SFTClusterContainer *CandCont,
			 const int *combination );

#endif
