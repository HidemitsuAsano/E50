/*
  ConfMan.hh

  2016/4  K.Shirotori
*/

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>

class DCGeomMan;

inline int FPhysProEMOn( int flag ){ return (flag | 0x01); }
inline int FPhysProDCYOn( int flag ){ return (flag | 0x02); }
inline int FPhysProHDOn( int flag ){ return (flag | 0x04); }

inline int FPhysProEMOff( int flag ){ return (flag & ~0x01); }
inline int FPhysProDCYOff( int flag ){ return (flag & ~0x02); }
inline int FPhysProHDOff( int flag ){ return (flag & ~0x04); }

inline bool GetFPhysProcEM( int flag ){ return (flag & 0x1); }
inline bool GetFPhysProcDCY( int flag ){ return (flag & 0x2); }
inline bool GetFPhysProcHD( int flag ){ return (flag & 0x4); }

inline int FStopNuOn( int flag ){ return (flag | 0x04); }
inline int FStopGamOn( int flag ){ return (flag | 0x08); }
inline int FStopEOn( int flag ){ return (flag | 0x10); }

inline int FStopNuOff( int flag ){ return (flag & ~0x04); }
inline int FStopGamOff( int flag ){ return (flag & ~0x08); }
inline int FStopEOff( int flag ){ return (flag & ~0x10); }

inline bool GetFStopNu( int flag ){ return (flag & 0x04); }
inline bool GetFStopGam( int flag ){ return (flag & 0x08); }
inline bool GetFStopE( int flag ){ return (flag & 0x10); }

class ConfMan
{
public:
  ConfMan( const std::string & filename );
  ~ConfMan();
private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );
public:
  static ConfMan *GetConfManager( void ) { return confManager_; }
  bool Initialize( void );
  void SetFileName( const std::string & filename ) { ConfFileName_=filename; }

  // Geometry
  DCGeomMan *GetDCGeomManager( void ) { return DCGeomManager_; }
  int GeomFlag( void ) const { return fGeom; }
  int DetectorType( void ) const { return fDetType_; }

  // Physics Process
  int PhysFlag( void ) const { return fPhysProc; }
  bool ExistEMProc( void ) const { return GetFPhysProcEM(fPhysProc); }
  bool ExistDCYProc( void ) const { return GetFPhysProcDCY(fPhysProc); }
  bool ExistHDProc( void ) const { return GetFPhysProcHD(fPhysProc); }

  // Stepping Action
  int StepFlag( void ) const { return fStepping; }
  bool DoesNuStop( void ) const { return GetFStopNu(fStepping); }
  bool DoesGamStop( void ) const { return GetFStopGam(fStepping); }
  bool DoesEStop( void ) const { return GetFStopE(fStepping); }

  // Primary Generator
  int ReactionMode( void ) const { return ReactionMode_; }
  double GetBeamX0( void ) const { return x0; }
  double GetBeamY0( void ) const { return y0; }
  double GetBeamZ0( void ) const { return z0; }
  double GetBeamU0( void ) const { return u0; }
  double GetBeamV0( void ) const { return v0; }
  double GetBeamP0( void ) const { return p0; }
  double GetBeamDX( void ) const { return dx; }
  double GetBeamDY( void ) const { return dy; }
  double GetBeamDZ( void ) const { return dz; }
  double GetBeamDU( void ) const { return du; }
  double GetBeamDV( void ) const { return dv; }
  double GetBeamDP( void ) const { return dp; }
  int GetPidBeam( void ) const { return beamPid; }
  bool MomBeam( void ) const { return MomBeam_; }

  // Sensitive Detector
  double GetTrEdep( void ) const { return trEdep_; }
  double GetCoEdep( void ) const { return coEdep_; }
  bool DetectEveryThing( void ) const { return detEVT_; }
  bool DetRootFill( void ) const { return detRootFill_; }

private:
  // ConfMan
  static ConfMan *confManager_;
  std::string ConfFileName_;

  // Geometry
  std::string DCGeomFileName_;
  DCGeomMan *DCGeomManager_;
  int fGeom;
  int fDetType_;

  // Physics Process
  int fPhysProc;

  // Steping Action
  int fStepping;

  // Primary Generation
  int ReactionMode_;

  // Beam parameters in mm, GeV, GeV/c and degree
  double x0,y0,z0,u0,v0,p0;
  double dx,dy,dz,du,dv,dp;
  bool MomBeam_;
  int beamPid;

  // Sensitive Detector
  double trEdep_, coEdep_;
  bool detEVT_;
  bool detRootFill_;

private:
  bool InitializeParameterFiles( void );
  void PrintParameters( void );
};

#endif
