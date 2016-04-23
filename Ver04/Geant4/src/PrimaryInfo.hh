/*
  PrimaryInfo.hh

  2012/5  K.Shirotori
*/

#ifndef PrimaryInfo_h
#define PrimaryInfo_h 1

struct PrimaryInfo 
{
  double x, y, z;
  double m1, m2;
  double p1, p2, p1cm, p2cm;
  double theta1, theta2, theta1cm, theta2cm;
  double phi1, phi2, phi1cm, phi2cm;
  double u1, u2, u1cm, u2cm;
  double v1, v2, v1cm, v2cm;
  double beta1, beta2;
  double pbeam, ubeam, vbeam;
  double p3, u3, v3;
  PrimaryInfo *Copy() const { return new PrimaryInfo(*this); }
};  

#endif
