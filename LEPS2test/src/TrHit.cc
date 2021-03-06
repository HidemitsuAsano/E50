//TrHits object has the vector of hit positions and TrLTrackHit in one segment (layer /fiber or ch)
//layer and ch information is stored when BSFTReco is called

#include <cmath>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "TrHit.hh"


TrHit::TrHit()
  : layer_(-1), fiber_(-1),adc_(-1),tot_(-1),localx_(-9999.),angle_(-100.)
{
  hitcounter_ = 0 ;
}

TrHit::TrHit( int layer, int ch )
  : layer_(layer), fiber_(ch),adc_(-1),tot_(-1),localx_(-9999.),angle_(-100.)
{
  hitcounter_ = 0 ;
}


TrHit::~TrHit()
{
 // clearRegisteredHits();
}



void TrHit::Print()
{
  std::cout << "TrHit::Print() " << std::endl;
  std::cout << "Layer:   " << layer_ << std::endl; 
  std::cout << "fiber:      " << fiber_ << std::endl; 
  std::cout << "hitcount " << hitcounter_ << std::endl;
  std::cout << "adc      " << adc_ << std::endl;
  std::cout << "tot      " << tot_ << std::endl;
  std::cout << "Local X  " << localx_ << std::endl;
  std::cout << "angle    " << angle_ << std::endl;
}




/*
void TrHit::SetPos( double pos )
{
  pos_.push_back(pos); 
  belongTrack_.push_back(false);
  dlRange_.push_back(true);
  hitcounter_++;
  hitID_.push_back(hitcounter_);
}*/

/*
void TrHit::clearRegisteredHits( void )
{
  int n=Cont_.size();
  for(int i=0; i<n; ++i)
    delete Cont_[i];
}
*/

//bool TrHit::CalcObservables( void )
//;{
//  static const std::string funcname="[TrHit::CalcObservables]";
  
//  ConfMan *confMan=ConfMan::GetConfManager();
//  if(!confMan){
//    std::cout << "Can not find ConfManager !! " << std::endl;
//    return false;
//  }
//
//  TrGeomMan *geomMan=confMan->GetTrGeomManager();
//  if(!geomMan){
//    std::cout << "Can not find GeomManager !! " << std::endl;
//    return false;
//  }
//   TrTdcCalibMan *calibMan=confMan->GetTrTdcCalibManager();
//   if(!calibMan) return false;
//   TrDriftParamMan *driftMan=confMan->GetTrDriftParamManager();
//   if(!driftMan) return false;

//  wpos_=geomMan->calcChPosition(layer_,ch_);
//  angle_=geomMan->GetTiltAngle(layer_);

//  std::cout << __FILE__ << " : " << __LINE__ << " layer: " << layer_ << " segment " <<ch_  <<  ": wpos_ "<< wpos_ << "angle " << angle_ << std::endl;
  
//  bool Status = true;
  //int nhitpos  = pos_.size();//nhitpos : number of hits 
  //for (int i=0; i<nhitpos; i++) {
    //     double ctime;
    //     if(!calibMan->GetTime( layer_, ch_, tdc_[i], ctime ))
    //       return false;
    
//     double dtime, dlength;
//     bool status=driftMan->calcDrift( layer_, ch_, ctime, dtime, dlength );
    
//     if (status == false) Status = status;
    
//     dt_.push_back(dtime);
    
    //Type 0, simple detector
    //if( confMan->AnaMode()==0 ){
    //  dl_.push_back(pos_[i]);
    //}//Type 1-3, realistic detector
   // else if( confMan->AnaMode()>=1 ){
    //  dl_.push_back(0.0);
   // }
    
//     if(layer_>=100){
//       if( dl_[i]>MinDLBc[layer_-100] && dl_[i]<MaxDLBc[layer_-100] )
// 	dlRange_[i]=true;
//     }
//     else{
//       if( dl_[i]>MinDLSdc[layer_] && dl_[i]<MaxDLSdc[layer_] )
// 	dlRange_[i]=true;
//     }
//   }
//  }
//
//  return Status;
//}
