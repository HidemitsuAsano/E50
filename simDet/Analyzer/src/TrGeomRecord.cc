/*
  TrGeomRecord.cc

  2012/5  K.Shirotori
*/

#include "ConfMan.hh"
#include "TrGeomRecord.hh"

#include <cmath>

const double Deg2Rad = acos(-1.)/180.;

//not used in the current code
void TrGeomRecord::calcVectors( void )
{
  //rotAngle1_ ,rotAngle2_: defined in a parameter file
  //tiltAngle : defined in a parameter file
  
  double ct1=cos(rotAngle1_*Deg2Rad), st1=sin(rotAngle1_*Deg2Rad);
  double ct2=cos(rotAngle2_*Deg2Rad), st2=sin(rotAngle2_*Deg2Rad);
  double ct0=cos(tiltAngle_*Deg2Rad), st0=sin(tiltAngle_*Deg2Rad);
  
  // rotation by Euler angle 
  // passive rotation
  //
  // |cos(psi) sin(psi)  0|
  // |-sin(psi) cos(psi) 0|   // rotAngle1 = psi (roation around z axis)
  // |0         0        1|
  // 
  // X
  //
  // |1           0           0| 
  // |0   cos(theta) sin(theta)|  //rotAngle2 = theta (rotation around x' axis)
  // |0  -sin(theta) cos(theta)|
  // 
  // X
  //
  // |cos(phi)  sin(phi) 0    |
  // |-sin(phi) cos(phi) 0    | //tiltAngle = phi (rotation around z'' axis)
  // |0         0        1    |
  //
  //    | dxds_  dyds_ dzds_ |
  // =  | dxdt_  dydt_ dzdt_ |
  //    | dxdu_  dydu_ dzdu_ |

  dxds_ =  ct0*ct2-st0*ct1*st2;
  dxdt_ = -st0*ct2-ct0*ct1*st2;
  dxdu_ =  st1*st2;

  dyds_ =  ct0*st2+st0*ct1*ct2;
  dydt_ = -st0*st2+ct0*ct1*ct2;
  dydu_ = -st1*ct2;

  dzds_ =  st0*st1;
  dzdt_ =  ct0*st1;
  dzdu_ =  ct1;

  dsdx_ =  ct0*ct2-st0*ct1*st2;
  dsdy_ =  ct0*st2+st0*ct1*ct2;
  dsdz_ =  st0*st1;

  dtdx_ = -st0*ct2-ct0*ct1*st2;
  dtdy_ = -st0*st2+ct0*ct1*ct2;
  dtdz_ =  ct0*st1;

  dudx_ =  st1*st2;
  dudy_ = -st1*ct2;
  dudz_ =  ct1;
}


//TODO : modify for type 2 and 3 detector
int TrGeomRecord::ChNumber( double pos ) const
{
  double dw=((pos-offset_)/dd_)+w0_;
  int iw=int(dw);
  if((dw-double(iw))>0.5)
    return iw+1;
  else
    return iw;
}

//Modified by H. Asano
//this function returns local x position of center of channel (fiber or a bundle of fibers)
//This calculation is for type 1 and 2 (round fiber)
double TrGeomRecord::ChPos (double ch) const
{ 
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  //(dd_ = dXdW) dd_ is the width of the fiber (or a bundle of fibers), 
  //dd_ = rSFT_Round_SegSpacing in GEANT4/r_DetectorSize.hh
  //offset_ is width/2.0 (mm) in type B and C detector
  //w0: offset for the first segment
  if( (confMan->AnaMode()) == 0 ){
    return dd_*(ch - w0_)+offset_; 
  }else if( (confMan->AnaMode()) > 0 ){
    //Type 1, 2 detector. 1 layer has 2 sublayers
    //dd_ : 
    //
    //return dd_*((int)((ch-1)/2.0+1) - w0_)+offset_*((int)(ch-1)%2);
    double pos = dd_*((int)((ch)/2.0+1) - w0_)+offset_*((int)ch%2);
     
    //if(Verbosity_){
    
    //{
    //  std::cout << __FILE__ << "  : " << __LINE__ << std::endl;
    //  std::cout << "AnaMode " << confMan->AnaMode() << std::endl;
    //  std::cout << "ch " << ch << std::endl;
    //  std::cout << "dd_ " << dd_ << std::endl;
    //  std::cout << "w0_ " << w0_ << std::endl;
    //  std::cout << "offset " << offset_ << std::endl;
    //  std::cout << "xl pos " << pos << std::endl;
    //}

    return pos;

  }else{
    std::cout << "Invalid AnaMode ! " << confMan->AnaMode() << std::endl;
    return -9999.;
  } 
}
