#ifndef GeomMan_h
#define GeomMan_h 1

#include <string>
#include <vector>
#include <map>

class GeomRecord;

class GeomMan 
{
private:
  GeomMan();
public:
  ~GeomMan();

public:
  void SetFileName( const char *filename ) { filename_=filename; }
  void SetFileName( const std::string &filename ) { filename_=filename; }

  bool Initialize( void );
  bool Initialize( const char *filename )
  { filename_=filename; return Initialize(); }
  bool Initialize( const std::string &filename )
  { filename_=filename; return Initialize(); }

  static GeomMan & GetInstance( void );
  double GetGlobalZ( int lnum ) const;
  double GetResolution( int lnum ) const;

  // Do not use this method except for special cases
  //void SetResolution( int lnum, double res ) const;
  //
  //void SetVertex( int lnum, ThreeVector vertex ) const;
  void SetChPattern( int pat ){ chpattern_ = pat;}

  double GetTiltAngle( int lnum ) const;
  double GetRotAngle1( int lnum ) const;
  double GetRotAngle2( int lnum ) const;
  double GetOffset( int lnum ) const;
  //const ThreeVector & GetGlobalPosition( int lnum ) const;
  //ThreeVector NormalVector( int lnum ) const;
  //ThreeVector UnitVector( int lnum ) const;
  const GeomRecord *GetRecord( int lnum ) const;

  //ThreeVector Local2GlobalPos( int lnum, const ThreeVector &in ) const;
  //ThreeVector Global2LocalPos( int lnum, const ThreeVector &in ) const;
  //ThreeVector Local2GlobalDir( int lnum, const ThreeVector &in ) const;
  //ThreeVector Global2LocalDir( int lnum, const ThreeVector &in ) const;

  double calcChPosition( int lnum, double ch ) const;
  int calcChNumber( int lnum, double position ) const;

  std::vector<int> GetDetectorIDList( void ) const;
  int GetDetectorId( const std::string &detName ) const;

  int GetSFTID(int lnum) const;
  int GetSFTLayer(int detectorID) const;
  
  //proto-type BFT's geoemetry

  int getXUV(int ch) const;
  int getsublayer(int ch) const ;
  int getlayer(int ch) const;
  int getbiglayer(int ch) const;
  int getfiber(int ch) const;
  int geteasiroc(int ch) const;
  int ChPattern(int ch) const; 

private:
  static GeomMan *geomMan_;
  std::string filename_;
  mutable std::map <int, GeomRecord *> geomRecord_;
  int chpattern_;
  
  void clearElements( void );
};


#endif
