/*
  Analysis.hh

  Base Class of Analysis

  2012/5  K.Shirotori
*/

#ifndef Analysis_h
#define Analysis_h 1

#include "G4String.hh"

#include <fstream>

// for gFile
//#include "TFile.h"

class G4Run;
class G4Event;
class PrimaryInfo;

class Analysis
{
public:
  Analysis( const G4String & histname );
  virtual ~Analysis();
private:
  Analysis( const Analysis & );
  Analysis & operator = ( const Analysis & );

public:
  virtual void BeginOfRun( const G4Run *aRun ) = 0;
  virtual void EndOfRun( const G4Run *aRun ) = 0;
  virtual void BeginOfEvent( const G4Event *anEvent ) = 0;
  virtual void EndOfEvent( const G4Event *anEvent ) = 0;
  virtual void PrimaryGeneration( const PrimaryInfo *info ) = 0;
  virtual void DefineHistograms( void ) = 0;
  //virtual void InitializeEvent(void) = 0;

  void SetFileName( const G4String &filename ) { filename_=filename; }
  G4bool GetTriggerStatus( void ) const { return fTriggered; }

  const G4String &GetFileName( void ) const { return filename_; }
  void SetActive( void ) { fActive_=true; }
  void SetInActive( void ) { fActive_=false; }

  void SaveFile( void ) const;
  void Terminate( void ) const;
  void ShowStatus( void ) const;

  void SetDataFile( const char *datafile );

protected:
  G4String filename_;
  G4bool fActive_;
  G4bool fTriggered;

  G4int trigNum;
  PrimaryInfo *pInfo_;
  std::ofstream DataFile_;

  virtual void PrintHitsInformation( const G4Event *anEvent,
				     std::ostream &ost ) const = 0;
};  

#endif
