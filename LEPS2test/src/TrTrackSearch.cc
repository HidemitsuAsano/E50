#include "TrTrackSearch.hh"
#include "TrParameters.hh"
//#include "TrLTrackHit.hh"
#include "SFTCluster.hh"
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
//const double MaxChisquare = 100000000.;//30
//const double MaxChisquareTr = 1000.;//30
const double MaxNumberOfClusters = 100.;//10.
const double MaxCombi = 1.0e6;

//const double MaxChisquareVXU = 1000.;
//const double ChisquareCutVXU = 1000.;
//int NhitGroup=0;
int Verbosity=0;
//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

//LocalTrackSearch 
//input : array of SFTClusterContainer , its index is layer number
//output :  vector of TrLocalTrack
//NumOfLayers : number of sft layers 12
//MinNumOfHits : numboer of required hits for locak tracking
//MaxChisquare : chi2/ndf cut of tracking. If this is <=0, no chi2/ndf cuts are applied 
int LocalTrackSearch(const  SFTClusterContainer *ClusterCont,
		      std::vector <TrLocalTrack *> &TrackCont, 
		      int NumOfLayers, unsigned int MinNumOfHits, double MaxChisquare )
{
  static const std::string funcname = "[LocalTrackSearch]";
  if(Verbosity>0){
    std::cout <<  __FILE__ << "  " << __LINE__ << "  " << funcname << std::endl;
  }

  const bool RequiredHitPattern[NumOfLayersSFT]=
  {
    true,//L0  X
    false,//L1  U
    false,//L2  V
    true,//L3  X
    false,//L4  U
    false,//L5  V
    false,//L6  U
    false,//L7  V
    true,//L8  X
    false,//L9  U
    false,//L10 V
    true,//L11 X
  };
  
  
  static int first=0;
  if(first==0){
    std::cout << std::endl;
    std::cout << __FILE__ << " L" << __LINE__ << " tracking parameters" << std::endl;
    std::cout << "# of layers " << NumOfLayers << std::endl;
    std::cout << "Minimum # of hits " << MinNumOfHits << std::endl;
    std::cout << "Max Chi2/ndf " << MaxChisquare << std::endl;
    std::cout << "Max # of clusters " <<  MaxNumberOfClusters << std::endl;
    std::cout << "Max # of Combination " << MaxCombi << std::endl;
    std::cout << "*** Required hit pattern*** " << std::endl;
    const char yes[2][4] = {"no","yes"};
    for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
      std::cout << "Layer" << ilr << "  "  << yes[RequiredHitPattern[ilr]] << std::endl;
    }
    std::cout << "********************" << std::endl;
    std::cout << std::endl;
    first++;
  }



  //vector of vector of SFTCluster 
  //1st row is layer 
  //2nd row is cluster in each layer
  //std::vector < std::vector <SFTCluster *> > CandCont;
  //std::vector <SFTCluster *>  CandCont[NumOfLayers];
  //CandCont.resize(NumOfLayers);
  
  //this way make too many candidates if there are multiple tracks (e.g. 2**12 layers) 
  int nCluster[NumOfLayers];
  for( int ilr=0; ilr<NumOfLayers; ++ilr ){ 
    nCluster[ilr]=(ClusterCont[ilr]).size();

    // If #Cluster>MaxNumerOfCluster,  error return
    if(nCluster[ilr]>MaxNumberOfClusters){
      std::cout << __FILE__ <<" : " << __LINE__ << " : "<< "MaxNumbeOfClusters exceed!! " << std::endl;
      std::cout << "layer: " << ilr << " MaxNumbefOfCluster : "  << MaxNumberOfClusters << " # of clusteres reconstructed : "
      << nCluster[ilr] << std::endl;
      for( int jlr=0; jlr<NumOfLayers; ++jlr )
        for_each( ClusterCont[jlr].begin(), ClusterCont[jlr].end(), DeleteObject() );
      return 0;
    } 
  }

#if 0
  std::cout << funcname << ": #Hits of each group" << std::endl;
  for( int i=0; i<NumOfLayers; ++i ) std::cout << std::setw(4) << nCluster[i];
  std::cout << std::endl;
  for( int i=0; i<NumOfLayers; ++i ){
    int n=ClusterCont[i].size();
    std::cout << "[" << std::setw(3) << i << "]: "
	      << std::setw(3) << n << " ";
    for( int j=0; j<n; ++j ){
      std::cout << ((TrLTrackHit *)ClusterCont[i][j]->GetHit(0))->GetWire() << " ";
    }
    std::cout << std::endl;
  }
#endif
  if(Verbosity>0){
    std::cout << __FILE__ << "   " << __LINE__ << "  making intex ...." << nCluster[0] << std::endl;
  }

  //make index of 2D vector (index of combination ,layer) . if the hit is associated to the track, the value of the component is -1 , if not 0
  std::vector < std::vector <int> > 
    CombiIndex = makeindex( NumOfLayers, 
			    MinNumOfHits, 
			    NumOfLayers, 
			    &(nCluster[0]) );
  int nnCombi=CombiIndex.size();
   
/*
  std::cout << __FILE__ << "   " << __LINE__ << " ===> " << nnCombi << " combinations will be checked.." 
    << std::endl;
  for(int ilr=0;ilr<NumOfLayers;ilr++){ 
    std::cout << "layer " << ilr << "  " <<  nCluster[ilr] << std::endl;
  }
  for(int icom = 0 ;icom < nnCombi ; icom++){
    int nlayer = CombiIndex[icom].size();
    for(int ilr=0;ilr<nlayer;ilr++){
      std::cout << CombiIndex[icom][ilr] << std::endl;
    }
  }
 */

#if 0
  std::cout << " ===> " << nnCombi << " combinations will be checked.." 
	    << std::endl;
#endif

  if( nnCombi>MaxCombi )  return 0;
  for( int i=0; i<nnCombi; ++i ){
    TrLocalTrack *track = MakeTrack( ClusterCont, &((CombiIndex[i])[0]) );
    if( !track ){std::cout << __FILE__ << "  " << __LINE__ << " no track " << std::endl; continue;}
    
    
    bool hitpattern[NumOfLayersSFT]={false};
    int nhit = track->GetNHit();
    //std::cout << __FILE__ << "  " << __LINE__ << "  " << nhit << std::endl;
    for(int ihit = 0 ; ihit< nhit; ihit++){
      SFTCluster *clsp = track->GetHit(ihit);
      if(clsp){
        int lnum = clsp->GetLayer();
        hitpattern[lnum] = true;
      }
      //std::cout << "ihit " << ihit << "layer " << clsp->GetLayer() << std::endl;
    }
    


    //for(int ilr=0;ilr<12;ilr++){
    //  std::cout << "layer: " << ilr << " Nhit " << track->GetHitOfLayerNumber(ilr) << std::endl;
    //}

    //check hit pattern
    bool isPatternOK=true;
    for( int ilr=0; ilr<NumOfLayersSFT; ilr++){
      if(RequiredHitPattern[ilr]){
        if(!hitpattern[ilr]){
          isPatternOK=false;
        }
      }
    }
    if(MaxChisquare<0) MaxChisquare = 99999999999.;

    if( track->GetNHit()>=MinNumOfHits && 
	track->DoFit() &&
	track->GetChiSquare()<MaxChisquare  &&
     isPatternOK
    ){
      TrackCont.push_back(track);
      //double chisqr = track->GetChiSquare();
    }
    else{
      /*
      *std::cout << "No tracks available" << std::endl;
      std::cout << "Nhit " << track->GetNHit() << std::endl;
      std::cout << "Fit? " << track->DoFit()  << std::endl;
      std::cout << "chi2 " << track->GetChiSquare() << std::endl;
      */
      delete track;
    }
  }

  // Clear Assciated track pointer
  int nbefore=TrackCont.size();
  for( int i=0; i<nbefore; ++i ){
    TrLocalTrack *tp=TrackCont[i];
    int nhit=tp->GetNHit(); 
    for( int j=0; j<nhit; ++j ) {
      tp->GetHit(j)->SetAssociatedLocalTrack(-1);//-1 is default value 
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
    for( int j=0; j<nhit; ++j ) tp->GetHit(j)->SetAssociatedLocalTrack(1);//temporality set 1
    for( int i2=TrackCont.size()-1; i2>i; --i2 ){
      TrLocalTrack *tp2=TrackCont[i2];
      int nhit2=tp2->GetNHit(), flag=0;
      for( int j=0; j<nhit2; ++j )
	if( tp2->GetHit(j)->GetAssociatedLocalTrack() ) ++flag;
      if(flag){
	delete tp2;
	TrackCont.erase(TrackCont.begin()+i2);
      }
    }      
  }
  
  //Set Track Id and  after deleting duplicated track above and
  //
  {
    int nn=TrackCont.size();
    for(int i=0; i<nn; ++i ){
      TrLocalTrack *tp=TrackCont[i];
      tp->SetTrackId(i);
      int nhit=tp->GetNHit();
      for( int j=0; j<nhit; ++j ){
	int lnum = tp->GetHit(j)->GetLayer();
  double zz = TrGeomMan::GetInstance().GetGlobalZ( lnum );
  tp->GetHit(j)->SetProjectedPosition(tp->GetX(zz), tp->GetY(zz));
  tp->GetHit(j)->SetAssociatedLocalTrack(i);
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

//  for( int i=0; i<NumOfLayers; ++i )
//    for_each( ClusterCont[i].begin(), ClusterCont[i].end(), DeleteObject() );
  
  return TrackCont.size();
}



//ndim = Number of layers
//ninimumHit : 8 is required so far
//index -> 12C12 + 12C11 + 12C10 + 12C9 + 12C8
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


//input :1. pointer of the  SFTClusterContainer (vector of cluster in each layer)
//       2. the number of combination for the track fit
//output: the pointer of the TrLocalTrack
//TrLocalTrack *MakeTrack(  const std::vector < std::vector <SFTCluster *> > &CandCont,
TrLocalTrack *MakeTrack(const SFTClusterContainer *CandCont,
			  const int *combination )
{
  static const std::string funcname = "[MakeTrack]";

  TrLocalTrack *tp=new TrLocalTrack();

  if(!tp){
    std::cerr << funcname << ": new fail" << std::endl;
    return 0;
  }    
  
  
  int nAssociatecluster[NumOfLayersSFT];
  for(int ilr=0;ilr<NumOfLayersSFT;ilr++){
    nAssociatecluster[ilr] = CandCont[ilr].size();
  }
  //int nAssociateCluster=CandCont->size();
  //std::cout << __FILE__ << "line " << __LINE__ << " number of associated cluster 0  " << nAssociateCluster[0] << std::endl;
  for( int ilr=0; ilr<NumOfLayersSFT; ++ilr ){
    int m=combination[ilr];
    SFTCluster *cluster=NULL;
    if(m>=0) cluster=CandCont[ilr][m];

///if 0
//    std::cout << funcname << ":" << std::setw(3)
//	      << iclus << std::setw(3) << m  << " " << std::endl;
	//      << CandCont[iclus][m] << std::endl; 
//#endif

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


    if(cluster){
      tp->AddHit( cluster );
    }
  }//iclus

  return tp;
}




