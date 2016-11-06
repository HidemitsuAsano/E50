/*
  ConfMan.cc

  2016/1  K.Shirotori
*/

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "PidUtil.hh"

#include <iostream>
#include <iomanip>
#include <sstream>

const std::string defDCGeomFile="DCgeom.param";
const std::string defFieldMapFile="FieldMap.dat";
const std::string defJAMEventGeneFile="JAMGeneData.dat";
const std::string defEvtGenDecayFile="param/EvtGenDecay.dat";
const std::string defEvtGenPDLFile="param/EvtGenPDL.dat";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename), DCGeomFileName_(defDCGeomFile),
    DCGeomManager_(0), TargSizX(-10.), TargSizY(-10.),TargSizZ(-10.),
    fGeom(0), fWCHole_(0), fDetType_(0), fIntImage_(0),
    matIdWorld(0), matIdDcGas1(0), matIdDcGas2(0), 
    matIdMagGap(0), matIdPMT(0), matIdPMTBox(0), matIdTarget(0), 
    matIdSMFAbsorber(0), matIdSP0Converter(0),
    fField_(false), FieldMapName_(defFieldMapFile), FieldScaleFactor(1.0), 
    fPhysProc(0), fStepping(0), ReactionMode_(0), 
    x0(0.0), y0(0.0), z0(0.0), u0(0.0), v0(0.0), p0(-1.0), 
    dx(0.0), dy(0.0), dz(0.0), du(0.0), dv(0.0), dp(0.0),
    GeneBeam_(false), MomBeam_(false),
    theta1(0.0), theta2(180.0), RecMass_(-1.), RecWidth_(-1.), CMAngle_(-1),
    beamPid(2), fData1_(false), fData2_(false), 
    JAMEventGeneName_(defJAMEventGeneFile), MaxJAMEvent_(0.0), DataCut_(-1),
    EvtGenDecayName_(defEvtGenDecayFile), 
    EvtGenPDLName_(defEvtGenPDLFile), 
    detEVT_(false), detRootFill_(false)
{
  static const std::string funcname="[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constring twice" << std::endl;
    exit(-1);
  }
  confManager_=this;
}

ConfMan::~ConfMan( )
{
  if(DCGeomManager_) delete DCGeomManager_;
  confManager_=0;
}

const G4int BufSize=300;

bool ConfMan::Initialize( void )
{
  static const std::string funcname="[ConfMan::Initialize]";

  char buf[BufSize], buf1[BufSize];
  G4double val1, val2, val3;
  G4int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){
      // Geometry
      if( sscanf(buf,"DCGEO: %s",buf1)==1 ){
	DCGeomFileName_=buf1;
      }
      else if( sscanf(buf,"TGTSIZ: %lf %lf %lf",&val1,&val2,&val3)==3 ){
	TargSizX=val1; 
	TargSizY=val2; 
	TargSizZ=val3;
      }
      else if( sscanf(buf,"TGTM: %d", &id)==1 ){
	if(id) fGeom=FGeomTargOn(fGeom);
	else   fGeom=FGeomTargOff(fGeom);
      }
      else if( sscanf(buf,"LQTGT: %d",&id)==1 ){
	if(id) fGeom=FGeomLiqTargOn(fGeom);
	else   fGeom=FGeomLiqTargOff(fGeom);
      }
      else if( sscanf(buf,"TGTLENGTH: %lf",&val1)==1 )
	TargSizZ=val1;
      else if( sscanf(buf,"WCHOLE: %d", &id )==1 ){
	if(id) fWCHole_=true;
	else   fWCHole_=false;
      }
      else if( sscanf(buf,"DETTYPE: %d", &id )==1 ){
	if(id) fDetType_=true;
	else   fDetType_=false;
      }
      else if( sscanf(buf,"INTIMAGE: %d", &id )==1 ){
	if(id) fIntImage_=true;
	else   fIntImage_=false;
      }
      // Materials
      else if( sscanf(buf,"MTWORLD: %d", &id)==1 )
	matIdWorld=id;
      else if( sscanf(buf,"MTDCGAS1: %d", &id)==1 )
	matIdDcGas1=id;
      else if( sscanf(buf,"MTDCGAS2: %d", &id)==1 )
	matIdDcGas2=id;
      else if( sscanf(buf,"MTMAGGAP: %d", &id)==1 )
	matIdMagGap=id;
      else if( sscanf(buf,"MTPMT: %d", &id)==1 )
	matIdPMT=id;
      else if( sscanf(buf,"MTPMTBOX: %d", &id)==1 )
	matIdPMTBox=id;
      else if( sscanf(buf,"MTTGT: %d", &id)==1 )
	matIdTarget=id;
      else if( sscanf(buf,"TGTMATER: %d", &id)==1 )
	matIdTarget=id;
      else if( sscanf(buf,"MTSMFABS: %d", &id)==1 )
	matIdSMFAbsorber=id;
      else if( sscanf(buf,"MTSP0CONV: %d", &id)==1 )
	matIdSP0Converter=id;
      // Field
      else if( sscanf(buf,"FLDMAP: %s",buf1)==1 ){
	FieldMapName_=buf1; fField_=true;
      }
      else if( sscanf(buf,"MAPSCALE: %lf", &val1)==1 )
	FieldScaleFactor=val1;
      // Physics Process
      else if( sscanf(buf,"EM: %d", &id )==1 ){
	if(id) fPhysProc=FPhysProEMOn(fPhysProc);
	else   fPhysProc=FPhysProEMOff(fPhysProc);
      }
      else if( sscanf(buf,"DECAY: %d", &id )==1 ){
	if(id) fPhysProc=FPhysProDCYOn(fPhysProc);
	else   fPhysProc=FPhysProDCYOff(fPhysProc);
      }
      else if( sscanf(buf,"HADRON: %d", &id )==1 ){
	if(id) fPhysProc=FPhysProHDOn(fPhysProc);
	else   fPhysProc=FPhysProHDOff(fPhysProc);
      }
      // Stepping Action
      else if( sscanf(buf,"SPECSTOP: %d", &id )==1 ){
	if(id) fStepping=FStopSPECOn(fStepping);
	else   fStepping=FStopSPECOff(fStepping);
      }
      else if( sscanf(buf,"SPECGSTOP: %d", &id )==1 ){
	if(id) fStepping=FStopSPECGamOn(fStepping);
	else   fStepping=FStopSPECGamOff(fStepping);
      }
      else if( sscanf(buf,"NUSTOP: %d", &id )==1 ){
	if(id) fStepping=FStopNuOn(fStepping);
	else   fStepping=FStopNuOff(fStepping);
      }
      else if( sscanf(buf,"GSTOP: %d", &id )==1 ){
	if(id) fStepping=FStopGamOn(fStepping);
	else   fStepping=FStopGamOff(fStepping);
      }
      else if( sscanf(buf,"ESTOP: %d", &id )==1 ){
	if(id) fStepping=FStopEOn(fStepping);
	else   fStepping=FStopEOff(fStepping);
      }
      // Primary Generation
      else if( sscanf(buf,"REACTION: %d", &id )==1 )
	ReactionMode_=id;
      else if( sscanf(buf,"BEAMX: %lf %lf",&val1, &val2 )==2 ){
	x0=val1; 
	dx=val2;
      }
      else if( sscanf(buf,"BEAMY: %lf %lf",&val1, &val2 )==2 ){
	y0=val1; 
	dy=val2;
      }
      else if( sscanf(buf,"BEAMZ: %lf %lf",&val1, &val2 )==2 ){
	z0=val1; 
	dz=val2;
      }
      else if( sscanf(buf,"BEAMU: %lf %lf",&val1, &val2 )==2 ){
	u0=val1; 
	du=val2;
      }
      else if( sscanf(buf,"BEAMV: %lf %lf",&val1, &val2 )==2 ){
	v0=val1; 
	dv=val2;
      }
      else if( sscanf(buf,"BEAMP: %lf %lf",&val1, &val2 )==2 ){
	p0=val1; 
	dp=val2;
      }
      else if( sscanf(buf,"GENEBEAM: %d", &id )==1 ){
	if(id) GeneBeam_=true;
	else   GeneBeam_=false;
      }
      else if( sscanf(buf,"MOMBEAM: %d", &id )==1 ){
	if(id) MomBeam_=true;
	else   MomBeam_=false;
      }
      else if( sscanf(buf,"ANGRANGE: %lf %lf",&val1, &val2 )==2 ){
	theta1=val1; 
	theta2=val2;
      }
      else if( sscanf(buf,"RMASS: %lf",&val1 )==1 )
	RecMass_=val1;
      else if( sscanf(buf,"RMWIDTH: %lf",&val1 )==1 )
	RecWidth_=val1;
      else if( sscanf(buf,"CMANGLE: %lf",&val1 )==1 )
	CMAngle_=val1;
      else if( sscanf(buf,"BPART: %d",&id )==1 )
	beamPid=id;
      else if( sscanf(buf,"MAXJAMEVENT: %lf",&val1 )==1 )
	MaxJAMEvent_=val1;
      else if( sscanf(buf,"JAMDATA: %s",buf1)==1 ){
	JAMEventGeneName_=buf1; 
	fData1_=true;
      }
      else if( sscanf(buf,"DATACUT: %d", &id )==1 ){
	if(id) DataCut_=true;
	else   DataCut_=false;
      }
      else if( sscanf(buf,"EVTGENDATA1: %s",buf1)==1 ){
	EvtGenDecayName_=buf1; 
	fEvtData1_=true;
      }
      else if( sscanf(buf,"EVTGENDATA2: %s",buf1)==1 ){
	EvtGenPDLName_=buf1; 
	fEvtData2_=true;
      }
      // Sensitive Detector
      else if( sscanf(buf,"DETECTALL: %d",&id )==1 ){
	if(id) detEVT_=true;
	else   detEVT_=false;
      }
      else if( sscanf(buf,"DETROOTFILL: %d",&id )==1 ){
	if(id) detRootFill_=true;
	else   detRootFill_=false;
      }
      else {
	std::cout << funcname << ": un-recognized record\n"
		  << buf << std::endl;
      }
    } 
  }

  fclose(fp);
  PrintParameters();
  InitializeParameterFiles();

  return true;
}


bool ConfMan::InitializeParameterFiles( void )
{
  DCGeomManager_ = & DCGeomMan::GetInstance();
  DCGeomManager_->Initialize(DCGeomFileName_);
  return true;
}

void ConfMan::PrintParameters( void )
{
  std::cout << "-----------------" << ConfFileName_ << "------" << std::endl;
  std::cout << "************** Geometry **************" << std::endl;
  std::cout << "DC Geom. Param.: " << DCGeomFileName_ << std::endl;
  std::cout << "HyperBall-J: " << GetFGeomHBJ(fGeom) << std::endl;
  std::cout << "Beamline: " << GetFGeomBL(fGeom) << std::endl;
  std::cout << "Target: " << GetFGeomTarg(fGeom) << std::endl;
  std::cout << "Target Size [mm]: " << TargSizX << " " << TargSizY
	    << " " << TargSizZ << std::endl;
  std::cout << "************* Materials *************" << std::endl;
  std::cout << "World: " << matIdWorld << std::endl;
  std::cout << "DCGas: " << matIdDcGas1 << " " << matIdDcGas2 << std::endl;
  std::cout << "MagGap: " << matIdMagGap << std::endl;
  std::cout << "PMT & PMTBox: " << matIdPMT << " " << matIdPMTBox << std::endl;
  std::cout << "SMF: " << matIdSMFAbsorber << std::endl;
  std::cout << "SP0: " << matIdSP0Converter << std::endl;
  std::cout << "Target: " << matIdTarget << std::endl;
  std::cout << "************ Field Map *************" << std::endl;
  std::cout << "Map File: " << FieldMapName_ << std::endl;
  std::cout << "Scale Factor: " << FieldScaleFactor << std::endl;
  std::cout << "************ Physics Process *******" << std::endl;
  std::cout << "EM:     " << GetFPhysProcEM(fPhysProc) << std::endl; 
  std::cout << "DECAY:  " << GetFPhysProcDCY(fPhysProc) << std::endl; 
  std::cout << "HADRON: " << GetFPhysProcHD(fPhysProc) << std::endl; 
  std::cout << "************ Stepping Action *******" << std::endl;
  std::cout << "SPEC Stop;           " << GetFStopSPEC(fStepping) << std::endl;
  std::cout << "SEPC gamma-ray Stop; " << GetFStopSPECGam(fStepping) << std::endl;
  std::cout << "Neutrino Stop;      " << GetFStopNu(fStepping) << std::endl;
  std::cout << "Gamma Stop;         " << GetFStopGam(fStepping) << std::endl;
  std::cout << "Electron Stop;      " << GetFStopE(fStepping) << std::endl;
  std::cout << "************ Primary Generation ****" << std::endl;
  std::cout << "Reaction Mode: " << ReactionMode_ << std::endl;
  std::cout << "Beam Particle: " << beamPid << " " 
	    << PIDParticleName(beamPid) << std::endl;
  std::cout << "Beam Param. X [mm]:    " << x0 << " " << dx << std::endl;
  std::cout << "Beam Param. Y [mm]:    " << y0 << " " << dy << std::endl;
  std::cout << "Beam Param. Z [mm]:    " << z0 << " " << dz << std::endl;
  std::cout << "Beam Param. U:         " << u0 << " " << du << std::endl;
  std::cout << "Beam Param. V:         " << v0 << " " << dv << std::endl;
  std::cout << "Beam Param. P [GeV/c]: " << p0 << " " << dp << std::endl;
  std::cout << "Ang. Range [degree]:   " << theta1 << " " 
	    << theta2  << std::endl;
  std::cout << "Recoil Particle Mass [GeV/c2]; " << RecMass_ << std::endl;
  std::cout << "Recoil Particle Width [MeV]; " << RecWidth_ << std::endl;
  std::cout << "JAM Generator File: " << JAMEventGeneName_ << std::endl;
  std::cout << "******** Sensitive Detector ****" << std::endl;
  std::cout << "Detect Everything: " << detEVT_ << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
}
