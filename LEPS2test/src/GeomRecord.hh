#ifndef GeomRecord_h
#define GeomRecord_h 1


#include <string>
#include <functional>

#include <TMath.h>
#include <Math/Vector3D.h>

class GeomRecord
{
public:
  GeomRecord( int id, const char *name,
                double x, double y, double z, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), offset_(ofs)
  { calcVectors(); }
    //perhaps, "length_" is confusing. This is just global z position defined in the parameter file
  GeomRecord( int id, const std::string &name,
                double x, double y, double z, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(x,y,z), tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), offset_(ofs)
  { calcVectors(); }
  
  
  GeomRecord( int id, const char *name,
                const ROOT::Math::XYZVector pos, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), offset_(ofs)
  { calcVectors(); }

  GeomRecord( int id, const std::string &name,
                const ROOT::Math::XYZVector pos, double ta,
                double ra1, double ra2, double length, double resol,
		double w0, double dd, double ofs )
    : id_(id), name_(name), pos_(pos),  tiltAngle_(ta),
      rotAngle1_(ra1), rotAngle2_(ra2),
      length_(length), resol_(resol), w0_(w0), dd_(dd), offset_(ofs)
  { calcVectors(); }
  
  ~GeomRecord() {}
  GeomRecord( const GeomRecord & );
  GeomRecord & operator=( const GeomRecord );

public:
  //const ThreeVector & Position( void ) const { return pos_; }
  //ThreeVector NormalVector( void ) const 
  //{ return ThreeVector( dxdu_, dydu_, dzdu_ ); }
  //ThreeVector UnitVector( void ) const
  //{ return ThreeVector( dxds_, dyds_, dzds_ ); }

  double dsdx( void ) const { return dsdx_; }
  double dsdy( void ) const { return dsdy_; }
  double dsdz( void ) const { return dsdz_; }
  double dtdx( void ) const { return dtdx_; }
  double dtdy( void ) const { return dtdy_; }
  double dtdz( void ) const { return dtdz_; }
  double dudx( void ) const { return dudx_; }
  double dudy( void ) const { return dudy_; }
  double dudz( void ) const { return dudz_; }

  double dxds( void ) const { return dxds_; }
  double dxdt( void ) const { return dxdt_; }
  double dxdu( void ) const { return dxdu_; }
  double dyds( void ) const { return dyds_; }
  double dydt( void ) const { return dydt_; }
  double dydu( void ) const { return dydu_; }
  double dzds( void ) const { return dzds_; }
  double dzdt( void ) const { return dzdt_; }
  double dzdu( void ) const { return dzdu_; }

  double ChPos( double ch ) const;

  //{ return dd_*(-ch + w0_)+offset_; }
  int ChNumber( double pos ) const; 

private:
  void calcVectors( void );

private:
  int id_;
  std::string name_;
  //ThreeVector pos_;
  ROOT::Math::XYZVector pos_;
  double tiltAngle_, rotAngle1_, rotAngle2_;
  double length_;
  double resol_;
  double w0_, dd_, offset_;

  double dxds_, dxdt_, dxdu_;
  double dyds_, dydt_, dydu_;
  double dzds_, dzdt_, dzdu_;

  double dsdx_, dsdy_, dsdz_;
  double dtdx_, dtdy_, dtdz_;
  double dudx_, dudy_, dudz_;
 
  friend class GeomMan;
  //friend class GeomRecordComp;

  //int Verbosity_;
};

/*
struct GeomRecordComp 
  : public std::binary_function <GeomRecord *, GeomRecord *, bool> 
{
  bool operator()( const GeomRecord * const p1,
		   const GeomRecord * const p2 ) const
  { return p1->id_ < p2->id_; }
};
*/

#endif
