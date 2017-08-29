#include <iostream>
#include <iomanip>

#include "SFTRawHit.hh"

SFTRawHit::SFTRawHit(int ch )
    : LayerId_(-1), ChId_(ch),type_(-1),
      //PosX_(0.0), PosY_(0.0), DL_(0.0),
      HitID_(-1),
      AdcHigh_(-1),
      otrAdcHigh_(false),
      AdcLow_(-1),
      otrAdcLow_(false),
      TdcLeading_(-1),
      TdcTrailing_(-1)
{


};


void SFTRawHit::Print()
{
  std::cout << std::endl;
  std::cout << "SFTRawHit::Print() " << std::endl;
  std::cout << "Layer ID  " << LayerId_ << std::endl;
  std::cout << "Ch ID     " << ChId_ << std::endl; 
  std::cout << "type      " << type_ << std::endl;
  std::cout << "fiber Id  " << fiberId_ << std::endl;
  std::cout << "Hit Id    " << HitID_ << std::endl;
  std::cout << "AdcHigh   " << AdcHigh_ << std::endl;
  std::cout << "otrAdcHigh" << otrAdcHigh_ << std::endl;
  std::cout << "AdcLow    " << AdcLow_ << std::endl;
  std::cout << "TDC lead  " << TdcLeading_ << std::endl;
  std::cout << "TDC trail " << TdcTrailing_ << std::endl;
}
