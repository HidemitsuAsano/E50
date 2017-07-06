//author: Hidemitsu Asano
//email: hidemitsu.asano@riken.jp
//date: June 23th, 2017


//how to Unpack NIMEasiroc raw data////////////////////////////////////
//1 word = 4 bytes = 32 bits and those 32bits are ordered as BigEndian.
//Then,
//1. run getBigEndian32() for 1 word 
//2. 32bits data should be shifted by Decode32bitWord();
//
//
//
//
//
//
//////////////////////////////////////////////////////////////////////


#ifndef RawData_h
#define RawData_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "DetectorID.hh"
#include <vector>

class PrimInfo;
class HodoRawHit;
class SFTRawHit;

typedef std::vector<PrimInfo*>   PrimInfoContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<SFTRawHit*>  SFTRawHitContainer;

class RawData
{

private:
  PrimInfoContainer PrimHC;

  HodoRHitContainer T0RHC;

  SFTRawHitContainer SFTRawHitCont;//vector of data object class for a single hit

public:
  RawData();
  ~RawData();

  void clearAll();
  bool Decode( std::ifstream & );
  unsigned int getBigEndian32(const char* b);
  unsigned int Decode32bitWord(unsigned int word32bit);
  bool isAdcHg(unsigned int data);
  bool isAdcLg(unsigned int data);
  bool isTdcLeading(unsigned int data);
  bool isTdcTrailing(unsigned int data);
  bool isScaler(unsigned int data);

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

  bool AddRawHit( SFTRawHitContainer &cont, 
		  int Ch, int hitID, 
       int AdcHigh, bool otrAdcHigh, 
       int AdcLow, bool otrAdcLow,
       int TdcLeading, int TdcTrailing);
  
  bool AddHodoRHit( HodoRHitContainer& cont,
		    int DetId, int Layer, int Seg,
		    double Time, double Edep, 
		    double Path, double Mom,
		    double PosX, double PosY, 
		    int Pid, double Beta );

public:

  const HodoRHitContainer& GetT0RHC() const;

  const SFTRawHitContainer & GetSFTRawHitContainer() const;
   
  //const SFTRawHit & GetSFTRawHit(int layer) const;
};

#endif

