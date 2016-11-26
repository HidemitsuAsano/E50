/*
  T0Hit.cc

  2012/5  K.Shirotori
*/

#include "T0Hit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <iomanip>

T0Hit::T0Hit()
  : layerID_(-1), segID_(-1), time_(-100.0), edep_(0.0),
    pos_(-9999.,-9999.,-9999.),
    mom_(-9999.,-9999.,-9999.),
    lmom_(-9999.,-9999.,-9999.),
    trackNo_(-1), fSignal_(false),
    path_(-9999.),beta_(-9999.),mass_(-9999.),
    Verbosity(0)
{}

G4Allocator<T0Hit> T0HitAllocator;

void T0Hit::Draw() const
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

void T0Hit::Print() const
{
  G4cout << "T0Hit " << std::setw(2) << layerID_ << "-" 
	 << std::setw(2) << segID_ << " > Time: " 
	 << time_/ns << " [ns] Edep.: " << edep_/MeV
	 << " [MeV]" << G4endl;
}

