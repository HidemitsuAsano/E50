#ifndef TrTrackSearch_h
#define TrTrackSearch_h 1

#include "TrAnalyzer.hh"
#include "DetectorID.hh"

#include <vector>

class SFTCluster;
class TrLocalTrack;
class TrLTrackHit;

//Full combination by Linear fitting
int LocalTrackSearch(const  SFTClusterContainer  *HC,
		      std::vector <TrLocalTrack *> &TrackCont,
		      int NumOfLayers, unsigned int MinNumOfHits, double MaxChisquare);

std::vector< std::vector<int> > makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 ); 
std::vector< std::vector<int> > makeindex_below( int ndim_org, int maximumHit, int ndim, const int *index1 ); 


TrLocalTrack *MakeTrack(const SFTClusterContainer *CandCont,
			 const int *combination );

// Structure, used only by TrTrackSearch
// which stores linked cluster information for SFT local tracking
// Since this class is only used by TrTrackSearch, it is
// locally stored in this file (i.e. this is a sort of Private Class)
struct SFTLinkedCluster {
  SFTLinkedCluster(SFTCluster *c, double xp, double yp, double magbend, int lay, int id):
    sftcluster_(c),projx_(xp),projy_(yp),mag_bend_(magbend), layer_(lay), clusterid_(id)
   {}

  
  SFTCluster *sftcluster_;
  double projx_;// projected position in global-x axis by straight line 
  double projy_;// projected position in global-y axis by straight line 
  double mag_bend_;//magnetic bend added to the straight line projection
  int layer_;//
  int clusterid_;//
};



struct SFTLinkedClusterContainer {
  SFTLinkedClusterContainer( ):
  nhit_(0),chi2_(-9999.),trackID_(-1){
    for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
      ishit_[ilr] = false;
    }
    for(int i=0;i<3;i++){
      p_[i] = -9999.;
    }
  }
  
  int nhit_;
  double chi2_;
  int trackID_;
  bool ishit_[NumOfLayersSFT];
  double p_[3];//3 momentum
  
  int GetNhit(void) {return nhit_;}
  bool IsHit(int layer) {return ishit_[layer];}
  void SetHit(int layer) { ishit_[layer] = true; }
  std::vector <SFTLinkedCluster> vlink;
};



#endif
