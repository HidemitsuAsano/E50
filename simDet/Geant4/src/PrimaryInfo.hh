/*
  PrimaryInfo.hh

  2012/5  K.Shirotori
*/

#ifndef PrimaryInfo_h
#define PrimaryInfo_h 1

struct PrimaryInfo 
{
  double x, y, z;//X,Y,Z primary vertex position
  double m1, m2;//beam mass ? , always 0
  double p1, p2, p1cm, p2cm;//magnitude of momentum, always 0
  double theta1, theta2, theta1cm, theta2cm;//angle of momentum, always 0
  double phi1, phi2, phi1cm, phi2cm;// always 0
  double u1, u2, u1cm, u2cm;// always 0
  double v1, v2, v1cm, v2cm;// always 0
  double beta1, beta2;// always 0
  //pbeam : magnitude of beam mom. ubeam : beamMom.x()/beamMom.z() vbeam : beamMom.y()/beamMom.z() 
  double pbeam, ubeam, vbeam;
  double xangle, yangle;//beam angle in x and y direction
  double p3, u3, v3;// ??? , always 0
  PrimaryInfo *Copy() const { return new PrimaryInfo(*this); }
};  

#endif
