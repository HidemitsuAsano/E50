
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "RawData.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "DetectorID.hh"
#include "SFTRawHit.hh"
#include "TemplateLib.hh"
#include "ConfMan.hh"
#include "GeomMan.hh"

RawData::RawData():
  SFTRawHitCont(0) //vector of data object class for a single track in SFT
{}

RawData::~RawData()
{
  clearAll();
}


bool RawData::AddRawHit( SFTRawHitContainer& cont,
			 int Ch, int hitID, 
       int AdcHigh, bool otrAdcHigh, 
       int AdcLow, bool otrAdcLow,
       int TdcLeading, int TdcTrailing)
{
  static const std::string funcname = "[RawData::AddRawHit]";
 
  SFTRawHit *p=0;
  int nh=cont.size();
  for( int i=0; i<nh; ++i ){
    SFTRawHit *q=cont[i];
    if( q->ChId()==Ch ){
      p=q; break;
    }
  }
  if(!p){
    p = new SFTRawHit(Ch);
    p->SetAdcHigh(AdcHigh);
    p->SetotrAdcHigh(otrAdcHigh);
    p->SetAdcLow(AdcLow);
    p->SetotrAdcLow(otrAdcLow);
    p->SetTdcLeading(TdcLeading);
    p->SetTdcTrailing(TdcTrailing);
    if(p) cont.push_back(p);
  }
  if(p){
    p->SetHitID( hitID );
    
    return true;
  }else{
    std::cerr << funcname << ": new fail." << std::endl;
    return false;
  }
}

void RawData::clearAll()
{
  std::for_each( SFTRawHitCont.begin(), SFTRawHitCont.end(), DeleteObject());
  SFTRawHitCont.clear();

  return;
}

//Decoding data 
//called event by event in ProcessingNormal(std::ifstream &In)
bool RawData::Decode( std::ifstream &In )
{
  clearAll();
  //ConfMan *confMan = ConfMan::GetConfManager();
  const GeomMan & geomMan=GeomMan::GetInstance();
  
  //# of total ch. hard-coded. modify later
  for(int ich=0;ich<128;ich++){
    SFTRawHit *p = new SFTRawHit(ich);
    SFTRawHitCont.push_back(p);
  }
  
  //number of NIM-EASIROCs
  const int nmodule = 2;

  for(int imod=0;imod < nmodule ;imod++){
    char headerByte[4];
    In.read(headerByte, 4);
    unsigned int header32 = getBigEndian32(headerByte);
    unsigned int header = Decode32bitWord(header32);
    bool isHeader = ((header >> 27) & 0x01) == 0x01;
    if(!isHeader) {
      std::cerr << __FILE__ << " L." << __LINE__ << std::endl;
      std::cerr << "Header data Frame Error" << std::endl;
      fprintf(stderr, "    %08X\n", header32);
      //std::exit(1);
      return false;
    }
    size_t dataSize = header & 0x0fff;
    //cout << dataSize << endl;
    
    //scaler is always off,
    //maybe, not impletemented in Chikuma-san's FPGA ? 
    //(version ver.20150923) 
    //unsigned int scalerValues[69];
    char* dataBytes = new char[dataSize * 4];
    In.read(dataBytes, dataSize * 4);
    for(size_t i = 0; i < dataSize; ++i) {
      unsigned int data32 = getBigEndian32(dataBytes + 4 * i);
      unsigned int data = Decode32bitWord(data32); 
      if(!data){
        std::cout << "Invalid data " << std::hex << data32 << std::dec << std::endl;
        std::cout << "Datasize " << dataSize << std::endl;
        return false;
      }
      int ch = (data >> 13) & 0x3f;
      if(imod%2==1) ch += 64; 
      SFTRawHitCont[ch]->SetChId(ch);
      int type = geomMan.getXUV(ch);
      SFTRawHitCont[ch]->SetType(type);
      int layer = geomMan.getlayer(ch);
      SFTRawHitCont[ch]->SetLayer(layer);
      //int fiber = geomMan.getfiber(ch); 
      //float pos = fiber/2.0 - 6.0;
      //int easiroc = geomMan.geteasiroc(ch);
      if(isAdcHg(data)) {
        bool otr = ((data >> 12) & 0x01) != 0;
        int value = data & 0x0fff;
        SFTRawHitCont[ch]->SetAdcHigh(value);
        SFTRawHitCont[ch]->SetotrAdcHigh(otr);
      }else if(isAdcLg(data)) {
        bool otr = ((data >> 12) & 0x01) != 0;
        int value = data & 0x0fff;
        SFTRawHitCont[ch]->SetAdcLow(value);
        SFTRawHitCont[ch]->SetotrAdcLow(otr);
      }else if(isTdcLeading(data)) {
        int value = data & 0x0fff;
        SFTRawHitCont[ch]->SetTdcLeading(value);
      }else if(isTdcTrailing(data)) {
        int value = data & 0x0fff;
        SFTRawHitCont[ch]->SetTdcTrailing(value);
      }else if(isScaler(data)) {
        //int value = data & 0x3fff;
        /*
        if(ch == 68) {
          int scalerValuesArrayIndex = events % 100;
          memcpy(scalerValuesArray[scalerValuesArrayIndex], scalerValues,
              sizeof(scalerValues));
        }*/
      }else {
        int ch = (data >> 13) & 0x3f;
        int value = data & 0x0fff;
        std::cout << __FILE__  << " l." <<  __LINE__ << " invalid data " << std::endl;
        std::cout << "adchg:"  << (data & 0x00680000);
        std::cout << "adclg:"  << (data & 0x00680000);
        std::cout << "tdcl:"   << (data & 0x00601000);
        std::cout << "tdct:"   << (data & 0x00601000);
        std::cout << "scaler:" << (data & 0x00600000);
        std::cout << "data:" << data << std::endl; 
        std::cout << "ch:" << ch << " value:" << value << std::endl;
        std::cerr << "Unknown data type" << std::endl;
      }
    }//end of loop of dataSize 
    delete[] dataBytes;
  }//for imodule
  

  return true;
}



const SFTRawHitContainer & RawData::GetSFTRawHitContainer() const
{
  return SFTRawHitCont;
}


unsigned int RawData::getBigEndian32(const char* b = NULL){
    //std::cout << "size of b " << sizeof(b) << std::endl;
    return ((b[0] << 24) & 0xff000000) |
           ((b[1] << 16) & 0x00ff0000) |
           ((b[2] <<  8) & 0x0000ff00) |
           ((b[3] <<  0) & 0x000000ff);
}


unsigned int RawData::Decode32bitWord(unsigned int word32bit = 0)
{
  //check data format
  unsigned int frame = word32bit & 0x80808080;
  if(frame != 0x80000000){
    std::cerr << __FILE__ << " " << __FUNCTION__ << " Frame Error! " << std::endl;
    std::cerr << "32 bit word: " << std::hex << word32bit << std::dec << std::endl;
    return 0;
  }

  return ((word32bit & 0x7f000000) >> 3) | 
         ((word32bit & 0x007f0000) >> 2) |
         ((word32bit & 0x00007f00) >> 1) |
         ((word32bit & 0x0000007f) >> 0);
}


//ADC High Gain
bool RawData::isAdcHg(unsigned int data = 0 )
{
    return (data & 0x00680000) == 0x00000000;
}


//ADC Low Gain
bool RawData::isAdcLg(unsigned int data = 0 )
{
    return (data & 0x00680000) == 0x00080000;
}

//TDC Leading
bool RawData::isTdcLeading(unsigned int data = 0 )
{
    return (data & 0x00601000) == 0x00201000;
}

//TDC Trailing
bool RawData::isTdcTrailing(unsigned int data = 0)
{
    return (data & 0x00601000) == 0x00200000;
}

//scaler function is not impletented yet in firmware (?) 
bool RawData::isScaler(unsigned int data = 0)
{
    return (data & 0x00600000) == 0x00400000;
}

