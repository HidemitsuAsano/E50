/*
  IITrackerHit.cc

  2013/4  K.Shirotori
*/

#include "IITrackerHit.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include <iomanip>

IITrackerHit::IITrackerHit()
  : layerID_(0), segID_(0), time_(-100.0), edep_(0.0),
    trackNo_(0), fSignal_(false)
{}

G4Allocator<IITrackerHit> IITrackerHitAllocator;

void IITrackerHit::Draw() const
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

void IITrackerHit::Print() const
{
  G4cout << "IITrackerHit " << std::setw(2) << layerID_ << "-" 
	 << std::setw(2) << segID_ << " > Time: " 
	 << time_/ns << " [ns] Edep.: " << edep_/MeV
	 << " [MeV]" << G4endl;
}

