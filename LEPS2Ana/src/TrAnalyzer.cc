#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <set>
#include <cmath>

#include "TrAnalyzer.hh"
#include "TrHit.hh"
#include "SFTCluster.hh"
#include "TrLocalTrack.hh"
#include "RawData.hh"
#include "SFTRawHit.hh"

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

TrAnalyzer::TrAnalyzer():
  isTrHitsSorted_(false),
  sigmathreshold_(10)
{
}

TrAnalyzer::~TrAnalyzer()
{
  clearTracksSFTT();
//  clearSFTClusters();
  clearTrHits();
}


bool TrAnalyzer::DecodeSFTRawHits( RawData *rawData )
{
  const std::string funcname = "[TrAnalyzer::DecodeSFTRawHits]";

  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan){
    std::cout << "Can not find ConfManager !! " << std::endl;
    return false;
  }

  TrGeomMan *geomMan=confMan->GetTrGeomManager();
  if(!geomMan){
    std::cout << "Can not find GeomManager !! " << std::endl;
    return false;
  }

  clearTrHits();
  clearSFTClusters();
  clearTracksSFTT();

  //SFT

  //Type 1,2,3 realistic detector
  static int first=0; 
  if( confMan->AnaMode()>=1 ){
    if(first==0){
      std::cout << __FILE__  << "  " << __LINE__ << " realistic detector is chosen " << std::endl;
      std::cout << __FILE__  << "  " << __LINE__ << " threshold of hits "  << sigmathreshold_ << " sigma" << std::endl;
      first++;
    }
    const SFTRawHitContainer &cont = rawData->GetSFTRawHitContainer();
    int nhit=cont.size();
      //std::cout<< nh << std::endl;
      
    for( int i=0; i<nhit; ++i ){
      SFTRawHit *rhit=cont[i];
      int rlayerID = rhit->LayerId();
      int rchID  = rhit->ChId();
      int adcval = rhit->GetAdcHigh();
      double pedestalmean = confMan->GetPedestalMean(rchID);
      double pedestalsigma = confMan->GetPedestalSigma(rchID);
      double thre = pedestalmean + sigmathreshold_*pedestalsigma;
      if(adcval < thre) continue;
      int fiber = geomMan->getfiber(rchID);
      //need to store fiber ID to sort hits for clusterig
      TrHit *hit = new TrHit(rlayerID ,fiber);
      int nhitpos= rhit->GetSize();
      for( int j=0; j<nhitpos; ++j ){
        hit->SetPos( fiber ); // set a segment ID as a hit position ??? ->re-fill in the CalcObservables
                          // May.23 2016 added comment: the TrHit object does not have the number of hits at this moment.Here, the vector of hit position is filled by setting the segment ID
      }
      //if(!hit) continue; 
	
      //should be use this function, tempolary hard-coded here. 
      //double wpos = geomMan->calcChPosition(rlayerID,rchID);
  
      int type = rhit->GetType();
      double offset = -9999;
      if(type==0){//x 
        offset = -3.0;
      }else{
        offset = -2.5;
      }
      double wpos = fiber/2.0 + offset;

      double angle = geomMan->GetTiltAngle(rlayerID);
#if check1
      std::cout<< __FILE__ << "  "  << __LINE__ << ": Layer  " << rhit->LayerId() << "  segment " << rhit->ChId() << 
        "  local x " << wpos << std::endl;
#endif
      hit->SetChPosition(wpos);
      hit->SetTiltAngle(angle);
  //std::cout << __FILE__ << " : " << __LINE__ << " layer: " << rlayerID << " segment " <<rchID  <<  ": wpos "<< wpos << "angle " << angle << std::endl;
  //if(hit->CalcObservables()){//hit position for each hit (before clustering) is calculated in this function, calling TrGeomRecord
      SFTTrHitContainer_[rlayerID].push_back(hit);
  //	}else{
  //   std::cout << __FILE__ << "TrHit::CalcObservables fail " << __LINE__ << std::endl;
  //  }
      }//for nhit
  }//if Type 1,2,3
  

  return true;
}


//sort TrHits for the clustering
//input: SFTTrHitContainer_ before sorting
//output: SFTTrHitContainer_ after sorting
bool TrAnalyzer::SortTrHits()
{
  //test
  /*
  for( int layer=0; layer<NumOfLayersSFT; ++layer ){
    int nhit = SFTTrHitContainer_[layer].size();
    std::cout << "layer : size" << layer << " : " << nhit << std::endl;
    for(int ihit = 0;ihit<nhit; ihit++){
      TrHit *hit = SFTTrHitContainer_[layer][ihit];
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer " << layer << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer stored " << hit->GetLayer() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "ch stored " << hit->GetCh() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "pos stored " << hit->GetChPosition() << std::endl; 
    }
  }*/
  //sort row hits here by segment id
  for( int layer=0; layer<NumOfLayersSFT; ++layer ){
    std::sort(SFTTrHitContainer_[layer].begin(),SFTTrHitContainer_[layer].end(),TrHit::compareTrHitPredicate);//
  }
  
  /*
  //test
  for( int layer=0; layer<NumOfLayersSFT; ++layer ){
    int nhit = SFTTrHitContainer_[layer].size();
    std::cout << "layer : size" << layer << " : " << nhit << std::endl;
    for(int ihit = 0;ihit<nhit; ihit++){
      TrHit *hit = SFTTrHitContainer_[layer][ihit];
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer " << layer << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "layer stored " << hit->GetLayer() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "ch stored " << hit->GetCh() << std::endl; 
      std::cout << "file " << __FILE__ << " line:  " << __LINE__ << "  " << "pos stored " << hit->GetPos() << std::endl; 
    }
  }*/

  isTrHitsSorted_ = true;

  return true;

}


////////////////////////////////////////////
//SFT clustering
////////////////////////////////////////////
//
//input: SFTTrHitContainer_ (member variable (vector of TrHit) of TrAnalyzer )
//output:SFTClusterContainer_ (member variable (vector of SFTCluster) of TrAnalyzer)
//MakeHitCluster is implemented in TrTrackSearch.cc as a non-member function
int TrAnalyzer::SFTClustering( void )
{
  for( int ilr=0; ilr<NumOfLayersSFT; ++ilr ){
    MakeHitCluster(SFTTrHitContainer_[ilr],SFTClusterContainer_[ilr]);
  }

  int nfoundcluster = 0;//TODO : not implemented yet
  return nfoundcluster;
}


//do clustering layer by layer
//input: TrHitContainer  (object of associated hits in each layer)
//output: the vector of SFTCluster 
bool TrAnalyzer::MakeHitCluster( const TrHitContainer &trhitcontainer,
     SFTClusterContainer &Cont )
{  
  ConfMan *confMan=ConfMan::GetConfManager();
  TrGeomMan *geomMan=confMan->GetTrGeomManager();
  int nhit=trhitcontainer.size(); //number of raw hits in the layer
  //std::cout << __FILE__ << "  " << __LINE__ << " nhit " << nhit << std::endl;
  if(nhit == 0) return true;
  std::vector <int> vLinkSegment; //store the candidate of segment # for clustering
  std::vector <double> vLinkSegmentPos; //store the mean position of candidate of segment # for clustering
  //std::vector <double> vLinkSegmentADC; //store the sum of ADC of candidate of segment # for clustering
  int clusterID=0;
  int vlinklxsize = 1;
  int vlinklzsize = 1;
  for( int ihit=0; ihit<nhit; ++ihit ){
    TrHit *hit=trhitcontainer[ihit];
    if( hit ){
      //int multiplicity = hit->GetPosSize();
      int layer = hit->GetLayer(); 
      int segment= hit->GetCh();
      double lxpos= hit->GetChPosition();//local x position
      /*
         if(multiplicity>1){
      //   std::cout << __FILE__ << "  " << __LINE__ << " multiple hits on one segment!!: " << multiplicity << std::endl;
      //   std::cout << "layer: " << layer << " segment: " << segment << std::endl;  
      }*/

      unsigned int vlinksize = vLinkSegment.size();
      bool isclusteringOK = false;
      
      //std::cout << __FILE__ << " : " << __LINE__ << " : " << "ihit " << ihit << " nhit " << nhit <<
      //" layer: " << layer << " segment " << segment << 
      //" local x pos " << lxpos << std::endl;
         //" clusterID: " << clusterID << std::endl;
         //std::cout << " size " << vlinksize << std::endl;
         
      if(vlinksize == 0){
        vLinkSegment.push_back(segment);
        vLinkSegmentPos.push_back(lxpos);
      }else{
        int seglinkcandidate1 = vLinkSegment.at(vlinksize-1);
        //candidates for clustering are always located in segment-1 or segment-2 
        //, since TrHits are sorted by ascending order
        if( ((segment - seglinkcandidate1) == 1) 
            || ((segment - seglinkcandidate1) == 2)
          ){
          //this cluster can be bigger, add this hit for clustering
          vLinkSegment.push_back(segment);
          vLinkSegmentPos.push_back(lxpos);
          //increase lzcluster size
          if( (segment - seglinkcandidate1) == 1) vlinklzsize++;
          
          //if there more than 1 hit for clustering, check also 2nd last hit
          int seglinkcandidate2 = -9999;//2nd last hit
          if(vlinksize>1) seglinkcandidate2 = vLinkSegment.at(vlinksize-2);
          
          if( (segment - seglinkcandidate1) == 2) vlinklxsize++;
          else if( (segment -seglinkcandidate2) == 2) vlinklxsize++;
        }else{
          isclusteringOK = true;

          /*
             std::cout << __FILE__ << " : " << __LINE__ << " clustering OK" << std::endl;

             std::cout << "ihit " << ihit << " nhit " << nhit << 
             " layer: " << layer <<
             " ID: " << clusterID << std::endl;
             std::cout << " size " << vlinksize << std::endl;
             */
        }//if clustering OK
      }//if vlinksize >0 

      if(isclusteringOK || (ihit == nhit-1 )){ // if it is last hit in the layer, finish clustering
        unsigned int currentvlinksize = vLinkSegment.size();

        SFTCluster *hitcluster = new SFTCluster();
        hitcluster->SetClusterID(clusterID);
        hitcluster->SetClusterSize(currentvlinksize);
        hitcluster->SetClusterLzSize(vlinklzsize);
        float LocalxSize  = (float) vlinklxsize;
        float offset = geomMan->GetOffset(layer);
        if(vlinklzsize%2 == 0) LocalxSize += offset;//offset from parameter file
        hitcluster->SetClusterLxSize(LocalxSize);

        double calclxpos=0.0;
        for(unsigned int jhit =0 ; jhit<currentvlinksize ;jhit++){
          calclxpos += vLinkSegmentPos.at(jhit);
        }
        //mean position of local-x
        calclxpos = calclxpos/(double)currentvlinksize;
        hitcluster->SetLocalX(calclxpos);
        double clsglobalz = geomMan->GetGlobalZ(layer);//get global z position from geometry file
        hitcluster->SetGlobalZ(clsglobalz);
        double tiltangle = geomMan->GetTiltAngle(layer);//degree 
        hitcluster->SetTiltAngle(tiltangle);
        hitcluster->SetLayer(layer);
        Cont.push_back(hitcluster);
           
           /*
           std::cout << __FILE__ << " : " << __LINE__ << " : " <<
           "ihit " << ihit << " nhit " << nhit << 
           " layer: " << layer <<
           " ID: " << clusterID << std::endl;
           std::cout << " cluster size: " << currentvlinksize << " lx cluster size " << LocalxSize << " lz cluster size " << vlinklzsize   << std::endl;
           std::cout << "clustered segment ";
           for(unsigned int icls =0 ;icls<currentvlinksize;icls++){
           std::cout << vLinkSegment.at(icls) << "  " ;
           }
           std::cout << std::endl;
           std::cout << "local x position " << calclxpos << std::endl;
           */

        //when clustering is finished,
        //clear the link vector and push back the next cluster candidate
        clusterID++;
        vLinkSegment.clear();
        vLinkSegmentPos.clear();
        vlinklxsize = 1;
        vlinklzsize = 1;
        vLinkSegment.push_back(segment);
        vLinkSegmentPos.push_back(lxpos);
      }//ifclustering OK
      // }//if Type 1,2,3 detector
    }//if TrHit
  }//for i hit

  return true;
}

//////////////////////////////////////////////////////
//Local tracking
/////////////////////////////////////////////////////

//SFT
bool TrAnalyzer::TrackSearchSFTT( void )
{
//  int ntrack =
    LocalTrackSearch( (SFTClusterContainer_), TrackSFTTCol, 
		      NumOfLayersSFT, MinNumOfHitsSFT );
//  std::cout<< __FILE__ << "  " << __LINE__ << "ntrack= " << ntrack << std::endl;
  
  return true;
}



void TrAnalyzer::clearTrHits( void )
{
  for( int l=0; l<NumOfLayersSFT; ++l ){
    for_each( SFTTrHitContainer_[l].begin(),  SFTTrHitContainer_[l].end(),  DeleteObject() );
    SFTTrHitContainer_[l].clear();
  }
}

void TrAnalyzer::clearSFTClusters( void )
{
  for( int l=0; l<NumOfLayersSFT; ++l ){
    for_each( SFTClusterContainer_[l].begin(),  SFTClusterContainer_[l].end(),  DeleteObject() );
    SFTClusterContainer_[l].clear();
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
