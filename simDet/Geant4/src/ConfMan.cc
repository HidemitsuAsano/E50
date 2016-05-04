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

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename), DCGeomFileName_(defDCGeomFile),
    DCGeomManager_(0), fGeom(0), fDetType_(0),
    fPhysProc(0), fStepping(0), ReactionMode_(0), 
    x0(0.0), y0(0.0), z0(0.0), u0(0.0), v0(0.0), p0(-1.0), 
    dx(0.0), dy(0.0), dz(0.0), du(0.0), dv(0.0), dp(0.0),
    MomBeam_(false), beamPid(1), 
    trEdep_(0.0), coEdep_(0.0),
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
  G4double val1, val2;
  G4int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }
  std::cout << "Conf file Open " << ConfFileName_.c_str() << std::endl;

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){
      // Geometry
      if( sscanf(buf,"DCGEO: %s",buf1)==1 ){
	DCGeomFileName_=buf1;
      } 
      else if( sscanf(buf,"DETTYPE: %d", &id )==1 ){
	fDetType_=id;
      std::cout << "reading DETTYPE " << fDetType_ << std::endl;
      }
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
      else if( sscanf(buf,"REACTION: %d", &id )==1 ){
	ReactionMode_=id;
      }
      else if( sscanf(buf,"BPART: %d",&id )==1 ){
	beamPid=id;
      }
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
      else if( sscanf(buf,"MOMBEAM: %d", &id )==1 ){
	if(id) MomBeam_=true;
	else   MomBeam_=false;
      }
      // Sensitive Detector
      else if( sscanf(buf,"TREDEP: %lf",&val1 )==1 ){
	trEdep_=val1; 
      }
      else if( sscanf(buf,"COEDEP: %lf",&val1 )==1 ){
	coEdep_=val1; 
      }
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
  std::cout << "************ Physics Process *******" << std::endl;
  std::cout << "EM:     " << GetFPhysProcEM(fPhysProc) << std::endl; 
  std::cout << "DECAY:  " << GetFPhysProcDCY(fPhysProc) << std::endl; 
  std::cout << "HADRON: " << GetFPhysProcHD(fPhysProc) << std::endl; 
  std::cout << "************ Stepping Action *******" << std::endl;
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
  std::cout << "******** Sensitive Detector ****" << std::endl;
  std::cout << "Detect Everything: " << detEVT_ << std::endl;
  std::cout << "--------------------------------------------" << std::endl;
}
