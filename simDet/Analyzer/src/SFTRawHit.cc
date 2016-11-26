/*
 SFTRawHit.cc

 2016/6  H.Asano
*/

#include "SFTRawHit.hh"


SFTRawHit::SFTRawHit( int layer, int ch )
    : LayerId_(layer), ChId_(ch), 
      PosX_(0.0), PosY_(0.0), DL_(0.0),
      HitID_(-1)
{


};


