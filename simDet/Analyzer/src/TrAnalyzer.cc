/*
  TrAnalyzer.cc

  2016/2  K.Shirotori
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>

#include "TrAnalyzer.hh"
#include "TrHit.hh"
#include "TrHitCluster.hh"
#include "TrLocalTrack.hh"
#include "RawData.hh"
#include "TrRawHit.hh"

#include "TemplateLib.hh"
#include "TrTrackSearch.hh"
#include "ConfMan.hh"
#include "TrParameters.hh"
#include "TrGeomMan.hh"
#include "DetectorID.hh"

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const double MaxChiSqrBeamTrack = 10000.;
const double MaxChiSqrScatTrack = 10000.;

#define check1 0

TrAnalyzer::TrAnalyzer()
{
}

TrAnalyzer::~TrAnalyzer()
{
  clearTracksSFTT();
  clearTrHitClusters();
  clearTrHits();
}


bool TrAnalyzer::DecodeRawHits( RawData *rawData )
{
  const std::string funcname = "[TrAnalyzer::DecodeRawHits]";

  ConfMan *confMan = ConfMan::GetConfManager();

  clearTrHits();
  clearTracksSFTT();

  //SFT

  //Type 0
  if( confMan->AnaMode()==0 ){
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetSFTRawHitContainer(layer);
      int nhit=cont.size();
      //std::cout<< nh << std::endl;
      
      for( int i=0; i<nhit; ++i ){
	TrRawHit *rhit=cont[i];
	
	TrHit *hit=new TrHit( rhit->LayerId(), rhit->WireId() );
	int nhitpos= rhit->GetSize();
	for( int j=0; j<nhitpos; ++j ){
	  hit->SetPos( rhit->GetDL(j) );
	  
#if check1
	  std::cout<< rhit->LayerId() << " " << rhit->GetDL(j) << std::endl;
#endif
	}
	//if(!hit) continue; 
	
	if(hit->CalcObservables())
	  SFTTrHitContainer_[layer].push_back(hit);
	else
	  delete hit;
      }
    }
  }//if Type 0

  //Type A , B , C
  if( confMan->AnaMode()>=1 ){
    for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
      const TrRHitContainer &cont =rawData->GetSFTRawHitContainer(layer);
      int nhit=cont.size();
      //std::cout<< nh << std::endl;
      
      for( int i=0; i<nhit; ++i ){
	TrRawHit *rhit=cont[i];
	
	TrHit *hit=new TrHit( rhit->LayerId(), rhit->WireId() );
	int nhitpos= rhit->GetSize();
  for( int j=0; j<nhitpos; ++j ){
	  hit->SetPos( rhit->WireId() ); // set a segment ID as a hit position ??? ->re-fill in the CalcObservables
                                   // May.23 2016 added comment: the TrHit object does not have the number of hits at this moment.Here, the vector of hit position is filled by setting the segment ID
	  
#if check1
	  std::cout<< __FILE__ << "  "  << __LINE__ << ": " << << rhit->LayerId() << " " << rhit->WireId() << std::endl;
#endif
	}
	//if(!hit) continue; 
	
	if(hit->CalcObservables()){//hit position for each hit (before clustering) is calculated in this function, calling TrGeomRecord
	  SFTTrHitContainer_[layer].push_back(hit);
	}else{
    std::cout << __FILE__ << "TrHit::CalcObservables fail " << __LINE__ << std::endl;
  }
  //  delete hit; //why ???
      }//for nhit
    }//for ilayer
  }//if Type A, B ,C
  
  

  return true;
}

bool TrAnalyzer::SortTrHits()
{

  //sort row hits here by wire id
  for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
    std::sort(SFTTrHitContainer_[layer].begin(),SFTTrHitContainer_[layer].end(),TrHit::compareTrHitPredicate);//simple sort does not work
  }

  //test
  for( int layer=1; layer<=NumOfLayersSFT; ++layer ){
    int nhit = SFTTrHitContainer_[layer].size();
    std::cout << "layer : size" << layer << " : " << nhit << std::endl;
    for(int ihit = 0;ihit<nhit; ihit++){
      TrHit *hit = SFTTrHitContainer_[layer][ihit];
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer " << layer << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer stored " << hit->GetLayer() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "wire stored " << hit->GetWire() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "pos stored " << hit->GetPos() << std::endl; 
    }
  }

  return true;

}


////////////////////////////////////////////
//SFT clustering
////////////////////////////////////////////
//
//input: SFTTrHitContainer_ (member variable (vector of TrHit) if TrAnalyzer )
//output:SFTTrHitClusterContainer_ (member variable (vector of TrHitCluster) if TrAnalyzer)
//MakeHitCluster is implemented in TrTrackSearch.cc as a non-member function
int TrAnalyzer::SFTClustering( void )
{
  for( int ilr=0; ilr<=NumOfLayersSFT; ++ilr ){
    MakeHitCluster(SFTTrHitContainer_[ilr],SFTTrHitClusterContainer_[ilr]);
  }

  int nfoundcluster = 0;//SFTTrHitClusterContainer_.size();TODO:how to get ??
  return nfoundcluster;
}

//////////////////////////////////////////////////////
//Local tracking
/////////////////////////////////////////////////////

//SFT
bool TrAnalyzer::TrackSearchSFTT( void )
{
  int ntrack =
    LocalTrackSearch( &(SFTTrHitContainer_[1]), TrackSFTTCol, 
		      NumOfLayersSFT, MinNumOfHitsSFT );
//  std::cout<< __FILE__ << "  " << __LINE__ << "ntrack= " << ntrack << std::endl;
  
  return true;
}



void TrAnalyzer::clearTrHits( void )
{
  for( int l=0; l<=NumOfLayersSFT; ++l ){
    for_each( SFTTrHitContainer_[l].begin(),  SFTTrHitContainer_[l].end(),  DeleteObject() );
    SFTTrHitContainer_[l].clear();
  }
}

void TrAnalyzer::clearTrHitClusters( void )
{
  for( int l=0; l<=NumOfLayersSFT; ++l ){
    for_each( SFTTrHitClusterContainer_[l].begin(),  SFTTrHitClusterContainer_[l].end(),  DeleteObject() );
    SFTTrHitClusterContainer_[l].clear();
  }
}

void TrAnalyzer::clearTracksSFTT( void )
{
  for_each( TrackSFTTCol.begin(), TrackSFTTCol.end(), DeleteObject() );
  TrackSFTTCol.clear();
}


bool TrAnalyzer::ReCalcTrHits( bool applyRecursively )
{
  // for( int l=0; l<=NumOfLayersSTIn; ++l ){
  //   int n=STInHC[l].size();
  //   for( int i=0; i<n; ++i ){
  //     TrHit *hit=(STInHC[l])[i];
  //     if(hit) hit->ReCalc(applyRecursively);
  //   }
  // }
 
  return true;
}

// bool TrAnalyzer::ReCalcTrackSTIn( bool applyRecursively )
// {
//   int n=TrackSTInCol.size();
//   for( int i=0; i<n; ++i ){
//     TrLocalTrack *track=TrackSTInCol[i];
//     if( track ) track->ReCalc( applyRecursively );
//   }
//   return true;
// }


bool TrAnalyzer::ReCalcAll( void )
{
  // ReCalcTrHits();

  // ReCalcTrackSTIn();

  return true;
}
