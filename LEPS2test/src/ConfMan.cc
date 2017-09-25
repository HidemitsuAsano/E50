#include "ConfMan.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include "GeomMan.hh"

const std::string defDCGeomFile="DCgeom.param";

ConfMan * ConfMan::confManager_ = 0;

ConfMan::ConfMan( const std::string & filename, const std::string tablename)
  : ConfFileName_(filename),
    PedetableName_(tablename),
    anaMode_(0), 
    GeomFileName_(defDCGeomFile),
    GeomManager_(0),
    SFTResol_(0.0),
    T0Resol_(0.0)
    //,ChMapPattern_(0)
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
  if(GeomManager_) delete GeomManager_;
  confManager_=0;
}

const int BufSize=1024;

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
	GeomFileName_=buf1;
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
  GeomManager_ = & GeomMan::GetInstance();
  //GeomManager_->SetChPattern(ChMapPattern_);
  GeomManager_->Initialize(GeomFileName_);
  std::ifstream ftable;
  ftable.open(PedetableName_.c_str());
  if(!ftable){
    std::cout << "NOT FOUND : " <<  PedetableName_.c_str() << std::endl;
    return false;
  }


  while(ftable){
    int ch;
    double mean,sigma;
    double hitmean,hitmeanerr;
    double hitsigma,hitsigmaerr;
    ftable >> ch >> mean >> sigma >> hitmean >> hitmeanerr >>  hitsigma >> hitsigmaerr;
    pedestalmean_[ch] = mean;
    pedestalsigma_[ch] = sigma;
    hitmean_[ch] = hitmean;
    hitmeanerr_[ch] = hitmeanerr;
    sigma_[ch] = hitsigma;
    sigmaerr_[ch] = hitsigmaerr;
  }

  return true;
}
