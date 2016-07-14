/*
  SpecField.cc

  2012/5  K.Shirotori
*/

#include "SpecField.hh"
#include "SimpleFieldElement.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"

SpecField::SpecField( const std::string &FieldMapName, double scaleFactor )
  : fMap( FieldMapName.c_str(),scaleFactor )
{
  fMap.Initialize();
}

SpecField::~SpecField()
{}

void SpecField::GetFieldValue( const double Point[4], 
			       double *Bfield ) const
{
  double X[3];
  X[0]=Point[0]/cm; X[1]=Point[1]/cm; X[2]=Point[2]/cm;

  if( fMap.GetFieldValue( X, Bfield ) ){
    Bfield[0] *= tesla;
    Bfield[1] *= tesla;
    Bfield[2] *= tesla;
  }
  else{
    Bfield[0]=Bfield[1]=Bfield[2]=0.0;
  }

  G4ThreeVector gPos( Point[0], Point[1], Point[2] );
  G4ThreeVector B( 0., 0., 0. );
  FMIterator end=elemList_.end();
  for( FMIterator itr=elemList_.begin(); itr!=end; ++itr ){
    if( (*itr)->ExistMagneticField() )
      B += (*itr)->GetMagneticField( gPos );
  }

  Bfield[0] += B.x(); Bfield[1] += B.y(); Bfield[2] += B.z();

#if 0
  G4cout << "X=(" << X[0] << "," << X[1] << "," << X[2] << ") "
	 << "B=(" << Bfield[0]/tesla << "," << Bfield[1]/tesla
	 << "," << Bfield[2]/tesla << ")" << G4endl;
#endif
}


void SpecField::cleanupSimpleElementList( void )
{
  elemList_.clear();
}

void SpecField::AddSimpleElement( SimpleFieldElement *elem )
{
  elemList_.push_back( elem );
}
