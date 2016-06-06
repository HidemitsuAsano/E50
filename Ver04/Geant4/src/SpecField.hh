/*
  SpecField.hh

  2012/5  K.Shirotori
*/

#ifndef SpecField_h 
#define SpecField_h 1

#include "SpecFieldMap.hh"
#include "G4MagneticField.hh"

#include <vector>

class SimpleFieldElement;

class SpecField : public G4MagneticField
{
public:
  SpecField();
  explicit SpecField( const std::string &FieldMapName, 
		      double scaleFactor=1.0 );
  ~SpecField();

private:
  SpecField( const SpecField & );
  SpecField & operator = ( const SpecField & );

public:
  void GetFieldValue( const double Point[4], double *Bfield ) const;
  void cleanupSimpleElementList( void );
  void AddSimpleElement( SimpleFieldElement *elem );

private:
  SpecFieldMap fMap;

  typedef std::vector <SimpleFieldElement *> FMContainer;
  typedef std::vector <SimpleFieldElement *>
  ::const_iterator FMIterator;

  FMContainer elemList_;
};

#endif

