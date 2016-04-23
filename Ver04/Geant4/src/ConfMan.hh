/*
  ConfMan.hh

  2016/1  K.Shirotori
*/

#ifndef ConfMan_h
#define ConfMan_h 1

#include <string>

class DCGeomMan;

inline int FGeomHBJOn( int flag ){ return  (flag | 0x02); }
inline int FGeomBLOn( int  flag ){ return (flag | 0x04); }
inline int FGeomTargOn( int  flag ){ return (flag | 0x08 ); }
inline int FGeomLiqTargOn( int  flag ){ return (flag | 0x10 ); }

inline int FGeomHBJOff( int flag ){ return  (flag & ~0x02); }
inline int FGeomBLOff( int  flag ){ return (flag & ~0x04); }
inline int FGeomTargOff( int  flag ){ return (flag & ~0x08 ); }
inline int FGeomLiqTargOff( int  flag ){ return (flag & ~0x10 ); }

inline bool GetFGeomHBJ( int flag ){ return (flag & 0x02); }
inline bool GetFGeomBL( int flag ){ return (flag & 0x04); }
inline bool GetFGeomTarg( int flag ){ return (flag & 0x08); }
inline bool GetFGeomLiqTarg( int flag ){ return (flag & 0x10); }

inline int FPhysProEMOn( int flag ){ return (flag | 0x01); }
inline int FPhysProDCYOn( int flag ){ return (flag | 0x02); }
inline int FPhysProHDOn( int flag ){ return (flag | 0x04); }

inline int FPhysProEMOff( int flag ){ return (flag & ~0x01); }
inline int FPhysProDCYOff( int flag ){ return (flag & ~0x02); }
inline int FPhysProHDOff( int flag ){ return (flag & ~0x04); }

inline bool GetFPhysProcEM( int flag ){ return (flag & 0x1); }
inline bool GetFPhysProcDCY( int flag ){ return (flag & 0x2); }
inline bool GetFPhysProcHD( int flag ){ return (flag & 0x4); }

inline int FStopSPECOn( int flag ){ return (flag | 0x01); }
inline int FStopSPECGamOn( int flag ){ return (flag | 0x02); }
inline int FStopNuOn( int flag ){ return (flag | 0x04); }
inline int FStopGamOn( int flag ){ return (flag | 0x08); }
inline int FStopEOn( int flag ){ return (flag | 0x10); }

inline int FStopSPECOff( int flag ){ return (flag & ~0x01); }
inline int FStopSPECGamOff( int flag ){ return (flag & ~0x02); }
inline int FStopNuOff( int flag ){ return (flag & ~0x04); }
inline int FStopGamOff( int flag ){ return (flag & ~0x08); }
inline int FStopEOff( int flag ){ return (flag & ~0x10); }

inline bool GetFStopSPEC( int flag ){ return (flag & 0x01); }
inline bool GetFStopSPECGam( int flag ){ return (flag & 0x02); }
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
  double GetTargetSizeX( void ) const { return TargSizX; }
  double GetTargetSizeY( void ) const { return TargSizY; }
  double GetTargetSizeZ( void ) const { return TargSizZ; }
  double TgtLength( void ) const { return TargSizZ; }
  int GeomFlag( void ) const { return fGeom; }
  bool ExistTarget( void ) const { return GetFGeomTarg(fGeom); };
  bool ExistLiqTarget( void ) const { return GetFGeomLiqTarg(fGeom); };
  bool WCHoleFlag( void ) const { return fWCHole_; }
  bool DetectorType( void ) const { return fDetType_; }

  // Materials
  int GetMatIdWorld( void ) const { return matIdWorld; }
  int GetMatIdDcGas1( void ) const { return matIdDcGas1; }
  int GetMatIdDcGas2( void ) const { return matIdDcGas2; }
  int GetMatIdMagGap( void ) const { return matIdMagGap; }
  int GetMatIdPMT( void ) const { return matIdPMT; }
  int GetMatIdPMTBox( void ) const { return matIdPMTBox; }
  int GetMatIdTarget( void ) const { return matIdTarget; }
  int GetMatIdSMFAbsorber( void ) const { return matIdSMFAbsorber; }
  int GetMatIdSP0Converter( void ) const { return matIdSP0Converter; }

  // Field
  bool ExistField( void ) const { return fField_; }
  const std::string & FieldMapName( void ) const { return FieldMapName_; }
  double MapScale( void ) const { return FieldScaleFactor; }

  // Physics Process
  int PhysFlag( void ) const { return fPhysProc; }
  bool ExistEMProc( void ) const { return GetFPhysProcEM(fPhysProc); }
  bool ExistDCYProc( void ) const { return GetFPhysProcDCY(fPhysProc); }
  bool ExistHDProc( void ) const { return GetFPhysProcHD(fPhysProc); }

  // Stepping Action
  int StepFlag( void ) const { return fStepping; }
  bool DoesStopInSPEC( void ) const { return GetFStopSPEC(fStepping); }
  bool DoesGamStopInSPEC( void ) const { return GetFStopSPECGam(fStepping); }
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
  bool GeneBeam( void ) const { return GeneBeam_; }
  bool MomBeam( void ) const { return MomBeam_; }

  double GetScatAngTheta1( void ) const { return theta1; }
  double GetScatAngTheta2( void ) const { return theta2; }
  double GetRecoilMass( void ) const { return RecMass_; }
  double GetRecoilWidth( void ) const { return RecWidth_; }
  int GetCMAngle( void ) const { return CMAngle_; }

  bool ExistData1( void ) const { return fData1_; }
  bool ExistData2( void ) const { return fData2_; }
  const std::string & JAMEventGeneName( void ) const { return JAMEventGeneName_; }
  double GetMaxJAMEvent( void ) const { return MaxJAMEvent_; }
  int DataCut( void ) const { return DataCut_; }
  bool ExistEvtData1( void ) const { return fEvtData1_; }
  bool ExistEvtData2( void ) const { return fEvtData2_; }
  const std::string & EvtGenDecayName( void ) const { return EvtGenDecayName_; }
  const std::string & EvtGenPDLName( void ) const { return EvtGenPDLName_; }

  // Sensitive Detector
  bool DetectEveryThing( void ) const { return detEVT_; }
  bool DetRootFill( void ) const { return detRootFill_; }

private:
  // ConfMan
  static ConfMan *confManager_;
  std::string ConfFileName_;

  // Geometry
  std::string DCGeomFileName_;
  DCGeomMan *DCGeomManager_;
  double  TargSizX, TargSizY, TargSizZ;
  int fGeom;
  bool fWCHole_, fDetType_;

  // Materials
  int matIdWorld, matIdDcGas1, matIdDcGas2, matIdMagGap,
    matIdPMT, matIdPMTBox, matIdTarget, matIdSMFAbsorber,
    matIdSP0Converter;

  // Field
  bool fField_;
  std::string FieldMapName_;
  double FieldScaleFactor;  // including polarity

  // Physics Process
  int fPhysProc;

  // Steping Action
  int fStepping;

  // Primary Generation
  int ReactionMode_;

  // Beam parameters in mm, GeV, GeV/c and degree
  double x0,y0,z0,u0,v0,p0;
  double dx,dy,dz,du,dv,dp;
  double theta1, theta2;
  double RecMass_; // c.f. Yc mass
  double RecWidth_; // c.f. Yc width
  int CMAngle_;
  int beamPid;
  bool GeneBeam_, MomBeam_;

  //JAM & Decay data files
  bool fData1_, fData2_;
  std::string JAMEventGeneName_;
  double MaxJAMEvent_;
  int DataCut_;

  //EvtGen data files
  bool fEvtData1_, fEvtData2_;
  std::string EvtGenDecayName_;
  std::string EvtGenPDLName_;

  // Sensitive Detector
  bool detEVT_;
  bool detRootFill_;

private:
  bool InitializeParameterFiles( void );
  void PrintParameters( void );
};


// Material Number (ID)
//  0 not defined 
//
//  1 Vacuum
//  2 Air
//  3 He Gas
//  4 Ar Gas
//  5 Ar+C2H6(50:50) Gas
//  6 At+iso-C4H10 (80:20) Gas
//  7 P10 Gas
//  
// 11 Fe
// 12 SUS316L
// 13 Al
// 14 Pb
// 15 Brass
//
// Followings are assumed only target material
// 21 Liq-H2
// 22 Liq-D2
// 23 Lig-He
// 24 Li
// 25 Be
// 26 B10
// 27 B11
// 28 C
// 29 LiO
// 30 LiN
// 31 Water


#endif
