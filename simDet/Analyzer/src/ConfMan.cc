/*
  ConfMan.cc

  2016/4 K.Shirotori
*/

#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "TrGeomMan.hh"

const std::string defDCGeomFile="DCgeom.param";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename )
  : ConfFileName_(filename),
    anaMode_(0), 
    TrGeomFileName_(defDCGeomFile),
    TrGeomManager_(0),
    SFTResol_(0.0),
    T0Resol_(0.0)
{
  static const std::string funcname="[ConfMan::ConfMan]";
  if( confManager_ ){
    std::cerr << funcname << ": constring twice" << std::endl;
    std::exit(-1);
  }
  confManager_=this;
}
ConfMan::~ConfMan( )
{
  if(TrGeomManager_) delete TrGeomManager_;
  confManager_=0;
}

const int BufSize=300;

bool ConfMan::Initialize()
{
  static const std::string funcname="[ConfMan::Initialize]";

  char buf[BufSize], buf1[BufSize];
  double val1;// val2, val3;
  int id;

  FILE *fp;
  if((fp=fopen(ConfFileName_.c_str(),"r"))==0){
    std::cerr << funcname << ": file open fail" << std::endl;
    exit(-1);
  }

  while( fgets( buf, BufSize, fp ) != 0 ){
    if( buf[0]!='#' ){ 
      //Analysis mode
      if( sscanf(buf,"ANAMODE: %d", &id )==1 ){
	anaMode_=id;
      }
      // Geometry
      else if( sscanf(buf,"TRGEO: %s",buf1)==1 ){
	TrGeomFileName_=buf1;
      }
      //Resolution
      else if( sscanf(buf,"SFTRESOL: %lf", &val1)==1 )
	SFTResol_=val1;
      else if( sscanf(buf,"T0RESOL: %lf", &val1)==1 )
	T0Resol_=val1;

      else {
	std::cout << funcname << ": un-recognized record\n"
		  << buf << std::endl;
      }
    }
  }

  fclose(fp);
  InitializeParameterFiles();
  InitializeHistograms();

  return true;
}

bool ConfMan::InitializeParameterFiles( void )
{
  TrGeomManager_ = & TrGeomMan::GetInstance();
  TrGeomManager_->Initialize(TrGeomFileName_);

  return true;
}
