/*
  TrGeomRecord.cc

  2012/5  K.Shirotori
*/

#include "TrGeomRecord.hh"

#include <cmath>

const double Deg2Rad = acos(-1.)/180.;

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


  dxds_ =  ct0*ct2-st0*ct1*st2;//always 0 in Ver.4
  dxdt_ = -st0*ct2-ct0*ct1*st2;//always 0 in Ver.4
  dxdu_ =  st1*st2;//-1 in Ver.4 (not always?)

  dyds_ =  ct0*st2+st0*ct1*ct2;//ct0 in ver.4
  dydt_ = -st0*st2+ct0*ct1*ct2;//st0 in ver.4
  dydu_ = -st1*ct2;//always 0 in Ver.4

  dzds_ =  st0*st1;//st0 in ver.4
  dzdt_ =  ct0*st1;//
  dzdu_ =  ct1;
  

  //Inverse matrix 
  // | dsdx_ dtdx_ dudx_ |
  // | dsdy_ dtdy_ dudy_ |
  // | dsdz_ dtdz_ dudz_ |

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

int TrGeomRecord::WireNumber( double pos ) const
{
  double dw=((pos-ofs_)/dd_)+w0_;
  int iw=int(dw);
  if((dw-double(iw))>0.5)
    return iw+1;
  else
    return iw;
}
