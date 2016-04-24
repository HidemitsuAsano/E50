/*
  HodoAnalyzer.hh

  2012/5  K.Shirotori
*/

#ifndef HodoAnalyzer_h
#define HodoAnalyzer_h

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class HodoHit;
class HodoCluster;

typedef std::vector <HodoHit*> HodoHitContainer;
typedef std::vector <HodoCluster*> HodoClusterContainer;

class HodoAnalyzer
{
public:
  HodoAnalyzer();
  ~HodoAnalyzer();
private:
  HodoAnalyzer(const HodoAnalyzer &);
  HodoAnalyzer & operator =(const HodoAnalyzer &);

private:
  HodoHitContainer T0Cont;

  HodoClusterContainer T0ClCont;

public:
  bool DecodeRawHits(RawData* rawData);

  bool DecodeT0Hits(RawData* rawData);  

  int GetNHitsT0( void ) const { return T0Cont.size(); };

  inline HodoHit * GetHitT0( int i ) const;

  int GetNClustersT1( void ) const { return T0ClCont.size(); };

  inline HodoCluster * GetClusterT0( int i ) const;

  bool ReCalcT0Hits( bool applyRecursively=false );
  
  bool ReCalcT0Clusters( bool applyRecursively=false );

  bool ReCalcAll( void );

private:
  void clearT0Hits();

  static int MakeUpClusters( const HodoHitContainer & HitCont,
			     HodoClusterContainer & ClusterCont,
			     double maxTimeDif );
};

inline HodoCluster * HodoAnalyzer::GetClusterT0( int i ) const
{
  if( i>=0 && i<T0ClCont.size() )
    return T0ClCont[i];
  else
    return 0;
}

inline HodoHit * HodoAnalyzer::GetHitT0( int i ) const
{
  if( i>=0 && i<T0Cont.size() )
    return T0Cont[i];
  else
    return 0;
}

#endif
