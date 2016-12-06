/*
  SFTHit.cc

  2016/4  K.Shirotori
*/

#include "SFTHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <iomanip>

SFTHit::SFTHit()
  : pass_(-9999),layerID_(-1), segID_(-1), time_(-100.0), edep_(0.0),
    pos_(-9999.,-9999.,-9999.),
    mom_(-9999.,-9999.,-9999.),
    lmom_(-9999.,-9999.,-9999.),
    trackNo_(0), fSignal_(false),
    xl_(-9999.),yl_(-9999.),
    xg_(-9999.),yg_(-9999.),zg_(-9999.),
    path_(-9999.),beta_(-9999.),mass_(-9999.)
{}

G4Allocator<SFTHit> SFTHitAllocator;

void SFTHit::Draw() const
{
  G4VVisManager *pVisManager = G4VVisManager::GetConcreteInstance();
  if( pVisManager ){
    G4Circle circle( pos_ );
    circle.SetScreenSize( 3.0 );              // in pixels
    circle.SetFillStyle( G4Circle::filled );
    G4Colour colour( 0., 1., 0. );            // green
    G4VisAttributes attribs( colour );
    circle.SetVisAttributes( attribs );
    pVisManager->Draw( circle );
  }
}

void SFTHit::Print() const
{
  G4cout << "SFTHit " << std::setw(2) << layerID_ << "-" 
	 << std::setw(2) << segID_ << " > Time: " 
	 << time_/ns << " [ns] Edep.: " << edep_/MeV
	 << " [MeV]" << G4endl;
}

