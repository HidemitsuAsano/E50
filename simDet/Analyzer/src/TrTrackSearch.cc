/*
  TrTrackSearch.cc

  2012/5  K.Shirotori
*/

#include "TrTrackSearch.hh"
#include "TrParameters.hh"
#include "TrLTrackHit.hh"
#include "TrHitCluster.hh"
#include "TrLocalTrack.hh"
#include "TemplateLib.hh"
#include "DetectorID.hh"
#include "TrGeomMan.hh"

#include "ConfMan.hh"

#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <cmath>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);
const double MaxChisquare = 100.;//30
const double MaxChisquareTr = 100.;//30
const double MaxNumberOfClusters = 100.;//10.
const double MaxCombi = 1.0e6;

const double MaxChisquareVXU = 100.;
const double ChisquareCutVXU = 100.;
int NhitGroup=0;

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////
int LocalTrackSearch( const TrHitContainer * hitcontainer,
		      std::vector <TrLocalTrack *> &TrackCont, 
		      int NumOfLayers, unsigned int MinNumOfHits )
{
  static const std::string funcname = "[LocalTrackSearch]";
  
  //vector of vector of TrHitCluster 
  //1st row is layer 
  //2nd row is cluster in each layer
  //std::vector < std::vector <TrHitCluster *> > CandCont;
  std::vector <TrHitCluster *>  CandCont[NumOfLayers];
  //CandCont.resize(NumOfLayers);
  
  //doing clustering here
  //TODO : move this to another function
  for( int ilr=0; ilr<NumOfLayers; ++ilr ){
    MakeHitCluster( hitcontainer[ilr], CandCont[ilr] );
  }
  
  int nCombi[NumOfLayers];
  for( int ilr=0; ilr<NumOfLayers; ++ilr ){ 
    nCombi[ilr]=(CandCont[ilr]).size();

    // If #Cluster>MaxNumerOfCluster,  error return
    if(nCombi[ilr]>MaxNumberOfClusters){
      std::cout << __FILE__ <<" : " << __LINE__ << " : "<< "MaxNumbeOfClusters exceed!! " << std::endl;
      std::cout << "layer: " << ilr << " MaxNumbefOfCluster : "  << MaxNumberOfClusters << " # of clusteres reconstructed : "<< 
    nCombi[ilr] << std::endl;
      for( int jlr=0; jlr<NumOfLayers; ++jlr )
        for_each( CandCont[jlr].begin(), CandCont[jlr].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCombi[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=CandCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)CandCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif
  
  std::vector < std::vector <int> > 
    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCombi[0]) );
  int nnCombi=CombiIndex.size();

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif
  if( nnCombi>MaxCombi )  return 0;
  for( int i=0; i<nnCombi; ++i ){
    TrLocalTrack *track = MakeTrack( CandCont, &((CombiIndex[i])[0]) );
    if( !track ) continue;
    if( track->GetNHit()>=MinNumOfHits && 
	track->DoFit() &&
	track->GetChiSquare()<MaxChisquare ){
      TrackCont.push_back(track);
  //    double chisqr = track->GetChiSquare();
    }
    else{
      //      std::cout << "No tracks available" << std::endl;
      delete track;
    }
  }

  // Clear Flags
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nhit=tp->GetNHit(); 
    for( int j=0; j<nhit; ++j ) {
      //if(tp->GetHit(j)->showFlags()) 
	tp->GetHit(j)->clearFlags();
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": Before Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  partial_sort( TrackCont.begin(), TrackCont.end(), 
		TrackCont.end(), TrLTrackComp() );

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Sorting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  // Delete Duplicated Tracks
  for( int i=0; i<int(TrackCont.size()); ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nhit=tp->GetNHit();
    for( int j=0; j<nhit; ++j ) tp->GetHit(j)->setFlags();

    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nhit2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nhit2; ++j )
	if( tp2->GetHit(j)->showFlags() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }

  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      int nhit=tp->GetNHit();
      for( int j=0; j<nhit; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
	double zz = TrGeomMan::GetInstance().GetLocalZ( lnum );
	tp->GetHit(j)->SetCalPosition(tp->GetX(zz), tp->GetY(zz));
      }
    }
  }

#if 0
  {
    int nn=TrackCont.size();
    std::cout << funcname << ": After Deleting. #Tracks = " 
	      << nn << std::endl;

    for( int i=0; i<nn; ++i ){
      TrLocalTrack *track=TrackCont[i];
      std::cout << std::setw(3) << i << " #Hits="
		<< std::setw(2) << track->GetNHit() 
		<< " ChiSqr=" << track->GetChiSquare()
		<< std::endl;
    }
    std::cout << std::endl;

  }
#endif

  for( int i=0; i<NumOfLayers; ++i )
    for_each( CandCont[i].begin(), CandCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}




std::vector< std::vector<int> > 
makeindex( int ndim_org, int minimumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex( ndim_org, minimumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if (validHitNum >= minimumHit)
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
    //  int size1=index.size();
    }
  }

  return index;
}

std::vector< std::vector<int> > 
makeindex_below( int ndim_org, int maximumHit, int ndim, const int *index1 )
{
  if(ndim==1){
    std::vector< std::vector<int> > index2;
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem(1,i);
      index2.push_back(elem);
    }
    return index2;
  }

  std::vector< std::vector<int> > 
    index2=makeindex_below( ndim_org, maximumHit, ndim-1, index1+1 );
 
  std::vector< std::vector<int> > index;
  int n2=index2.size();
  for( int j=0; j<n2; ++j ){
    for( int i=-1; i<index1[0]; ++i ){
      std::vector <int> elem;
      int validHitNum=0;
      int n3=index2[j].size();
      elem.reserve(n3+1);
      elem.push_back(i);
      if (i != -1)
	validHitNum++;
      for( int k=0; k<n3; ++k ) {
        elem.push_back(index2[j][k]);
        if (index2[j][k] != -1)
          validHitNum++;
      }
      if (ndim==ndim_org) {
        if ((validHitNum <= maximumHit) && (validHitNum>0))
          index.push_back(elem);
      } else {
        index.push_back(elem);
      }
      //int size1=index.size();
    }
  }

  return index;
}


//make cluster layer by layer
//input: TrHitContainer 
//output: the vector of TrHitCluster 
bool MakeHitCluster( const TrHitContainer & trhitcontainer,
		     std::vector <TrHitCluster *> & Cont )
{  
//  ConfMan *confMan=ConfMan::GetConfManager();
//  std::cout << __FILE__ << std::endl;
  int nhit=trhitcontainer.size(); //number of raw hits in the layer
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
      /*
         if( confMan->AnaMode()==0 ){
         for (int im=0; im<multiplicity; im++) {
      //if( !(hit->rangecheck(m)) ) continue; //checking the range of drift length

      //double pos=hit->GetPos(im);	
      //double wp=hit->GetWirePosition();//local-x position
      //double dl=hit->GetDriftLength(m);

      //Cont.push_back( new TrHitCluster( new TrLTrackHit(hit,pos,im) ) );
      }
      }else if( confMan->AnaMode()>=1 ){//Type A, B,C detector */

      //int layer = hit->GetLayer(); 
      int segment= hit->GetWire();
      double lxpos=hit->GetWirePosition();//local-x position
      /*
         if(multiplicity>1){
      //   std::cout << __FILE__ << "  " << __LINE__ << " multiple hits on one segment!!: " << multiplicity << std::endl;
      //   std::cout << "layer: " << layer << " segment: " << segment << std::endl;  
      }*/

      unsigned int vlinksize = vLinkSegment.size();
      int seglinkcandidate = -9999;
      bool isclusteringOK = false;
      /*
         std::cout << __FILE__ << " : " << __LINE__ << " : " << "ihit " << ihit << " nhit " << nhit <<
         " layer: " << layer << " segment " << segment << 
         " clusterID: " << clusterID << std::endl;
         std::cout << " size " << vlinksize << std::endl;
         */
      if(vlinksize == 0){
        vLinkSegment.push_back(segment);
        vLinkSegmentPos.push_back(lxpos);
      }else{
        seglinkcandidate = vLinkSegment.at(vlinksize-1);
        //std::cout << "segment " << segment << " linkcandidate " << seglinkcandidate << std::endl;
        //candidates for clustering are always located in segment-1 or segment-2 
        //, since TrHits are sorted by ascending order
        if( ((segment - seglinkcandidate) == 1) 
            || ((segment - seglinkcandidate) == 2)
          ){
          //this cluster can be bigger
          //add this hit for clustering
          vLinkSegment.push_back(segment);
          vLinkSegmentPos.push_back(lxpos);
          //std::cout << "add segment " << segment << std::endl;
          //increase lzcluster size
          if( (segment - seglinkcandidate) == 1) vlinklzsize++;
          //increase lxcluster size
          if( (segment - seglinkcandidate) == 2) vlinklxsize++;
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

      if(isclusteringOK || (ihit == nhit-1 )){ // if it is last hit, finish clustering
        unsigned int currentvlinksize = vLinkSegment.size();

        TrHitCluster *hitcluster = new TrHitCluster();
        hitcluster->SetClusterID(clusterID);
        hitcluster->SetClusterSize(currentvlinksize);
        hitcluster->SetClusterLzSize(vlinklzsize);
        float LocalxSize  = (float) vlinklxsize;
        if(vlinklzsize%2 == 0) LocalxSize += 0.5;//this should be offset from param.file
        hitcluster->SetClusterLzSize(LocalxSize);

        double calclxpos=0;
        for(unsigned int jhit =0 ; jhit<currentvlinksize ;jhit++){
          calclxpos += vLinkSegmentPos.at(jhit);
        }
        //mean position of local-x
        calclxpos = calclxpos/(double)currentvlinksize;
        hitcluster->SetLocalX(calclxpos);


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
           std::cout << "local -x " << calclxpos << std::endl;
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
      // }//if Type A, B, C detector
    }//if TrHit
  }//for i hit

  return true;
}

//TrLocalTrack *MakeTrack(  const std::vector < std::vector <TrHitCluster *> > &CandCont,
TrLocalTrack *MakeTrack(  std::vector <TrHitCluster *> *CandCont,
			  const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  TrLocalTrack *tp=new TrLocalTrack();

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    

  //int n=CandCont.size();
  const int n = 12;//NumOfLayers;
  for( int i=0; i<n; ++i ){
    int m=combination[i];
    TrHitCluster *cluster=0;
    if(m>=0) cluster=CandCont[i][m];
#if 0
    std::cout << funcname << ":" << std::setw(3)
	      << i << std::setw(3) << m  << " "
	      << CandCont[i][m] << std::endl; 
#endif

#if 0
    /*****************************************************************/
    // Only SdcIn & BcOut Local Tracking !!                          
    /*****************************************************************/
    int mm = 0 ;
    //std::cout << "-------------------------------------" << std::endl;
    if(cluster){
     mm = cluster->NumberOfHits();
     if(mm>2) std::cout << "mm = " << mm << std::endl;

     double DL[mm] ;
     double DT[mm] ;
     int Layer[mm] ;

      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);

	DL[j] = hitp->GetDriftLength();
	DT[j] = hitp->GetDriftTime();
	Layer[j] = hitp->GetLayer();

	if(mm<2 || !((Layer[j]>0 && Layer[j]<11) || (Layer[j]>112 && Layer[j]<125))){
	  if((Layer[j]>=1 && Layer[j]<=4) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=5 && Layer[j]<=10) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=113 && Layer[j]<=118) && (DT[j]>-10 && DT[j]<40))
	    if(hitp) tp->AddHit( hitp );
	  if((Layer[j]>=119 && Layer[j]<=124) && (DT[j]>-10 && DT[j]<60))
	    if(hitp) tp->AddHit( hitp );
	}
      }
      //std::cout << "***********************************" << std::endl;
      if(mm==2 && ((Layer[0]>0 && Layer[0]<11) || (Layer[0]>112 && Layer[0]<125))){
	int Ok = 0 ;
	/*	
	std::cout << "DL[ " << Layer[0] << "] = " << DL[0] 
		  << " DL[ " << Layer[1] << "] = " << DL[1] << std::endl; 
	*/
	//std::cout << "***********************************" << std::endl;
	if(Layer[0] ==1 || Layer[0]==3 || Layer[0]==113 || Layer[0]==115 || Layer[0]==117){
	  if(((DL[0]+DL[1])>0.7) && ((DL[0]+DL[1])<2.3)) Ok = 1 ;
	}  
	if(Layer[0] ==5 || Layer[0]==7 || Layer[0]==9 || Layer[0]==119 || Layer[0]==121 || Layer[0]==123){
	  if(((DL[0]+DL[1])>1.5) && ((DL[0]+DL[1])<3.5)) Ok = 1 ;
	}
	if(Ok){
	  for(int j=0; j<mm; ++j ){
	    TrLTrackHit *hitp=cluster->GetHit(j);
	    if(hitp) tp->AddHit( hitp );
	  }
	}
      }
    }
#endif

/*
#if 1
    if(cluster){
     int mm=cluster->NumberOfHits();
      for(int j=0; j<mm; ++j ){
	TrLTrackHit *hitp=cluster->GetHit(j);
	  if(hitp) tp->AddHit( hitp );
      }
    }
#endif
*/
  }
  return tp;
}

