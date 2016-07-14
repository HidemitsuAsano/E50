/*
  ReadJAMData.hh

  2013/4  K.Shirotori
*/

#ifndef ReadJAMData_h
#define ReadJAMData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cstdio>

#include <vector>
#include "G4ThreeVector.hh"

class JAMData;

typedef std::vector<JAMData*> JAMDataContainer;

class ReadJAMData
{
public:
  ReadJAMData();
  ~ReadJAMData();

private:
  ReadJAMData(const ReadJAMData&);
  ReadJAMData& operator=(const ReadJAMData&);
  
private:
  JAMDataContainer Data;
  G4int evnum_;

public:
  void clearAll();
  bool Processing( std::istream &, G4int EventNum );

private:
  bool AddData( JAMDataContainer& cont,
		G4int EvNum, G4int Pid, 
		G4double MomX, G4double MomY, G4double MomZ );

  
public:
  const JAMDataContainer& GetJAMData() const;

};

#endif

