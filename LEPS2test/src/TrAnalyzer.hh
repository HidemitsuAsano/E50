#ifndef TrAnalyzer_h 
#define TrAnalyzer_h 1

//author: H.Asano
//email: hidemitsu.asano@riken.jp
//last update Aug.13th,2017
//class name: TrAnalyzer
//
//Main class to reconstruct event and called in AnaBFT.cc 
//This class works for 
//1. settting threshold of ADC and TDC. defining hits
//2. clustering
//3. calling several functions for tracking

#include "DetectorID.hh"
#include <vector>

class TrHit;
class SFTCluster;
class TrLocalTrack;
class RawData;

//vector of hit/cluster object in each layer 
typedef std::vector <TrHit *> TrHitContainer;
typedef std::vector <SFTCluster *> SFTClusterContainer;//

class TrAnalyzer
{
public:
  TrAnalyzer();
  ~TrAnalyzer();
private:
  TrAnalyzer( const TrAnalyzer & );
  TrAnalyzer & operator = ( const TrAnalyzer & );
  bool isTrHitsSorted_;
  double sigmathreshold_;
  bool isTOTcut_;
  double TOTcutlow_;
  double TOTcutup_;
  TrHitContainer SFTTrHitContainer_[NumOfLayersSFT];
  SFTClusterContainer SFTClusterContainer_[NumOfLayersSFT];

  std::vector <TrLocalTrack *> TrackSFTTCol;//vector of SFT track class 
  //chi2/ndf cut for tracking
  double MaxChisquare_;

public:
  bool MakeSFTRawHits( RawData *rawData );
  bool SortTrHits();
  void SetSigmaThreshold(double sigma) { sigmathreshold_ = sigma; };
  double GetSigmaThreshold() {return sigmathreshold_;};
  void UseTOTcut(bool yes=true){ isTOTcut_ = yes;};
  void SetTOTcut(double low,double up);
  inline const TrHitContainer & GetSFTTrHitContainer( int layer ) const;
  inline const SFTClusterContainer & GetSFTClusterContainer( int layer ) const;
   
  int SFTClustering(void);
  bool MakeHitCluster(const TrHitContainer &trhitcontainer,SFTClusterContainer &Cont );
  bool TrackSearchSFTT( void );

  int GetNtracksSFTT( void ) const  { return TrackSFTTCol.size(); }

  inline TrLocalTrack * GetTrackSFTT(unsigned int i ) const;

  bool ReCalcTrHits( bool applyRecursively=false ); 
  bool ReCalcTrackSFTT( bool applyRecursively=false ); 

  bool ReCalcAll( void );
  void SetMaxChisuquare( double cut) { MaxChisquare_ = cut;}

private:
  void clearTrHits( void );
  void clearSFTClusters( void );
  void clearTracksSFTT( void );

public:
  void resetTracksSFTT( void ) { clearTracksSFTT(); }
};

inline const TrHitContainer & TrAnalyzer::GetSFTTrHitContainer( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=-1;
  return SFTTrHitContainer_[layer];
}

inline const SFTClusterContainer & TrAnalyzer::GetSFTClusterContainer( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=-1;
  return SFTClusterContainer_[layer];
}

inline TrLocalTrack * TrAnalyzer::GetTrackSFTT( unsigned int i ) const
{
  if(  i<TrackSFTTCol.size() )
    return TrackSFTTCol[i];
  else
    return 0;
}

#endif 
