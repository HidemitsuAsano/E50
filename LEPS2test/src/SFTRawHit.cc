#include "SFTRawHit.hh"


SFTRawHit::SFTRawHit(int ch )
    : LayerId_(-1), ChId_(ch),type_(-1),
      PosX_(0.0), PosY_(0.0), DL_(0.0),
      HitID_(-1),
      AdcHigh_(-1),
      otrAdcHigh_(false),
      AdcLow_(-1),
      otrAdcLow_(false),
      TdcLeading_(-1),
      TdcTrailing_(-1)
{


};


