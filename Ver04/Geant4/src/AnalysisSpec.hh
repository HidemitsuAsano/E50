/*
  AnalysisSpec.hh

  2012/5  K.Shirotori
*/

#ifndef AnalysisSpec_h
#define AnalysisSpec_h 1

#include "Analysis.hh"

class AnalysisSpec : public Analysis
{
public:
  AnalysisSpec( const G4String & histname );
  ~AnalysisSpec();

public:
  virtual void BeginOfRun( const G4Run *aRun );
  virtual void EndOfRun( const G4Run *aRun );
  virtual void BeginOfEvent( const G4Event *anEvent );
  virtual void EndOfEvent( const G4Event *anEvent );
  virtual void PrimaryGeneration( const PrimaryInfo *info );
  virtual void InitializeEvent(void);
  virtual void DefineHistograms( void );

private:
  virtual void PrintHitsInformation( const G4Event *anEvent, 
				     std::ostream &ost ) const;
};

#endif
