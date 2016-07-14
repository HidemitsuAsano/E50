/*
  Analysis.cc

  Base Class of Analysis

  2012/5  K.Shirotori
*/

#include "Analysis.hh"
#include "PrimaryInfo.hh"
#include "G4ios.hh"
#include "TFile.h"

#include <iomanip>

#include <signal.h>

void SigHandler( G4int )
{
  if(gFile){
    gFile->Write();
    gFile->Close();
  }
}

Analysis::Analysis( const G4String & histname )
  : filename_(histname), fActive_(true), fTriggered(false),
    pInfo_(0), DataFile_(0)
{
  signal( SIGINT, SigHandler );
  signal( SIGABRT, SigHandler );
}

Analysis::~Analysis()
{
  SaveFile();
  if(pInfo_) delete pInfo_;
}


void Analysis::SaveFile( void ) const
{
  if( fActive_ )
    gFile->Write();
}

void Analysis::Terminate( void ) const
{
  G4cout << "[Analysis] Terminate()" << G4endl;   
  if( fActive_ ){
    gFile->Write();
    gFile->Close();
  }
}

void Analysis::ShowStatus( void ) const
{
  G4cout << "Analyzer Status\n"
         << "  File : " << filename_ << "\n"
         << "  Analysis : " << fActive_ << "\n"
         << "  Trigger : " << fTriggered << "\n"
         << "  # of triggered events : " << trigNum 
         << G4endl;
}

void Analysis::SetDataFile( const char *datafile )
{
  DataFile_.open( datafile );
}
