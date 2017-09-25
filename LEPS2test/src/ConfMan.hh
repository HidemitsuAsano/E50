#ifndef Confman_h
#define Confman_h 1

#include <string>

class VEvent;
class GeomMan;
class FieldMan;

class ConfMan
{

private:
  //Conf
  std::string ConfFileName_;
  std::string PedetableName_;
  static ConfMan *confManager_;

  //Analysis mode
  int anaMode_;

  //Geometry files  
  std::string GeomFileName_;
  GeomMan *GeomManager_;

  //Tracker position resolution
  double SFTResol_;

  //Counter time resolution
  double T0Resol_;
  
  //not used now
  //int ChMapPattern_;
  
  double pedestalmean_[128];
  double pedestalsigma_[128];
  double hitmean_[128];
  double hitmeanerr_[128];
  double sigma_[128];
  double sigmaerr_[128];

public:
  explicit ConfMan( const std::string & filename,const std::string table );
  ~ConfMan();

  static ConfMan *GetConfManager( void ) { return confManager_; }
  bool Initialize( void );
  void SetFileName( const std::string & filename ) { ConfFileName_=filename; }
  VEvent* EventAllocator();

  int AnaMode( void ) const { return anaMode_; }

  GeomMan *GetGeomManager( void ) { return GeomManager_; }
  //void SetChMap(int pat) { ChMapPattern_ = pat;}

  double GetSFTResol( void ) const { return SFTResol_; }

  double GetT0Resol( void ) const { return T0Resol_; }
  double GetPedestalMean(int ch) { return pedestalmean_[ch];}
  double GetPedestalSigma(int ch) { return pedestalsigma_[ch];}
  double GetHitMean(int ch) { return hitmean_[ch];}
  double GetHitMeanError(int ch) { return hitmeanerr_[ch];}
  double GetHitSigma(int ch) { return sigma_[ch];}
  double GetHitSigmaError(int ch) { return sigmaerr_[ch];}

private:
  ConfMan( const ConfMan & );
  ConfMan & operator = ( const ConfMan & );

public:
  bool InitializeHistograms();
  bool InitializeParameterFiles();

};

#endif
