/*
  TrAnalyzer.hh

  2016/2  K.Shirotori
*/

#ifndef TrAnalyzer_h 
#define TrAnalyzer_h 1

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

class TrHit;
class SFTCluster;
class TrLocalTrack;
class RawData;

//vector of hit/cluster object in each layer 
typedef std::vector <TrHit *> TrHitContainer;
typedef std::vector <SFTCluster *> SFTClusterContainer;//added by H.Asano

class TrAnalyzer
{
public:
  TrAnalyzer();
  ~TrAnalyzer();
private:
  TrAnalyzer( const TrAnalyzer & );
  TrAnalyzer & operator = ( const TrAnalyzer & );
  bool isTrHitsSorted_;
private:
  TrHitContainer SFTTrHitContainer_[NumOfLayersSFT];
  SFTClusterContainer SFTClusterContainer_[NumOfLayersSFT];

  std::vector <TrLocalTrack *> TrackSFTTCol;//vector of SFT track class 

public:
  bool DecodeSFTRawHits( RawData *rawData );
  bool SortTrHits();

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
  if( i>=0 && i<TrackSFTTCol.size() )
    return TrackSFTTCol[i];
  else
    return 0;
}

#endif 
