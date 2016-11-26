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
  double ct1=cos(rotAngle1_*Deg2Rad), st1=sin(rotAngle1_*Deg2Rad);
  double ct2=cos(rotAngle2_*Deg2Rad), st2=sin(rotAngle2_*Deg2Rad);
  double ct0=cos(tiltAngle_*Deg2Rad), st0=sin(tiltAngle_*Deg2Rad);

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
int TrGeomRecord::WireNumber( double pos ) const
{
  double dw=((pos-offset_)/dd_)+w0_;
  int iw=int(dw);
  if((dw-double(iw))>0.5)
    return iw+1;
  else
    return iw;
}

//Modified by H. Asano
//this function returns local position of center of fiber (or ch)
//The calculation is different between type 1,2 and 3 detector
//
double TrGeomRecord::WirePos (double ch) const
{ 
  ConfMan *confMan=ConfMan::GetConfManager();
  if(!confMan) return false;
  //(dd_ = dXdW) dd_ is probably the width of the fiber
  //offset_ is width/2.0 (mm) in type B and C detector
  //w0: offset for the first segment
  if( (confMan->AnaMode()) <= 1 ){
    return dd_*(ch - w0_)+offset_; 
  }else if( (confMan->AnaMode()) > 1 ){//Type B, C detector. 1 layer has 2 sublayer
    //return dd_*((int)((ch-1)/2.0+1) - w0_)+offset_*((int)(ch-1)%2);
    double pos = dd_*((int)((ch)/2.0+1) - w0_)+offset_*((int)ch%2);
    /*
    std::cout << "ch " << ch << std::endl;
    std::cout << "dd_ " << dd_ << std::endl;
    std::cout << "w0_ " << w0_ << std::endl;
    std::cout << "offset " << offset_ << std::endl;
    std::cout << "xl pos " << pos << std::endl;
    */
    return pos;
  }else{
    return -9999.;
  } 
}
