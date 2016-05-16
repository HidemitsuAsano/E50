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
class TrLocalTrack;
class RawData;

typedef std::vector <TrHit *> TrHitContainer;

class TrAnalyzer
{
public:
  TrAnalyzer();
  ~TrAnalyzer();
private:
  TrAnalyzer( const TrAnalyzer & );
  TrAnalyzer & operator = ( const TrAnalyzer & );

private:
  TrHitContainer SFTTHC[NumOfLayersSFT+1];

  std::vector <TrLocalTrack *> TrackSFTTCol;//array of SFT track class 

public:
  bool DecodeRawHits( RawData *rawData );

  inline const TrHitContainer & GetSFTTHC( int layer ) const;

  bool TrackSearchSFTT( void );

  int GetNtracksSFTT( void ) const  { return TrackSFTTCol.size(); }

  inline TrLocalTrack * GetTrackSFTT( int i ) const;

  bool ReCalcTrHits( bool applyRecursively=false ); 
  bool ReCalcTrackSFTT( bool applyRecursively=false ); 

  bool ReCalcAll( void );

private:
  void clearTrHits( void );
  void clearTracksSFTT( void );

public:
  void resetTracksSFTT( void ) { clearTracksSFTT(); }
};

inline const TrHitContainer & TrAnalyzer::GetSFTTHC( int layer ) const
{
  if( layer<0 || layer>NumOfLayersSFT ) layer=0;
  return SFTTHC[layer];
}

inline TrLocalTrack * TrAnalyzer::GetTrackSFTT( int i ) const
{
  if( i>=0 && i<TrackSFTTCol.size() )
    return TrackSFTTCol[i];
  else
    return 0;
}

#endif 
