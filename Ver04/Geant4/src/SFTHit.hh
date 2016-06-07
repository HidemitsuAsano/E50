/*
  SFTHit.hh

  2016/4  K.Shirotori
*/

#ifndef SFTHit_h
#define SFTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"

class SFTHit : public G4VHit
{
public:
  SFTHit();
  ~SFTHit() { }
private:
  SFTHit( const SFTHit & );
  SFTHit & operator = ( const SFTHit & );
public:
  int operator == ( const SFTHit & ) const { return 0; }

  inline void * operator new ( size_t size );
  inline void operator delete( void * aHit );

  void Draw() const;
  void Print() const;

private:
  G4int pass_;
  G4int layerID_; 
  G4int segID_;
  G4double time_;
  G4double edep_;
  G4ThreeVector pos_;
  G4ThreeVector mom_;
  G4ThreeVector lmom_;
  G4int trackNo_;
  G4bool fSignal_;
  G4double xl_, yl_;
  std::vector <G4String> HitPartName_;
  G4double path_;
  G4double beta_;
  G4double mass_;

public:
  void SetPass( void ) { pass_=1; }
  void SetLayerID( G4int id ) { layerID_=id; }
  void SetSegmentID( G4int id ) { segID_ = id; }
  void SetTime( G4double time ) { time_=time; }
  void SetEdep( G4double edep ) { edep_ = edep; }
  void AddEdep( G4double edep ) { edep_ += edep; }
  void SetPos( const G4ThreeVector &pos ) { pos_=pos; }
  void SetMom( const G4ThreeVector &mom ) { mom_=mom; }
  void SetLMom( const G4ThreeVector &mom ) { lmom_=mom; }
  void SetTrackNo( G4int no ) { trackNo_=no; }
  void SetTrueSignal() { fSignal_=true; }
  void SetFalseSignal() { fSignal_=false; }
  void SetLocalPos( G4double x, G4double y ) { xl_=x; yl_=y; }
  void SetHitParticleName( const G4String & name ) 
  { HitPartName_.push_back( name ); }
  void SetPathLength( G4double path ) { path_=path; }
  void SetBeta( G4double beta ) { beta_=beta; }
  void SetHitParticleMass( G4double mass ) { mass_=mass; }

  G4int GetPass() const { return pass_; }
  G4int GetLayerID( void ) const { return layerID_; }
  G4int GetSegmentID() const { return segID_; }
  G4double GetTime( void ) const { return time_; }
  G4double GetEdep() const { return edep_; }
  G4ThreeVector GetPos( void ) const { return pos_; }
  G4ThreeVector GetMom( void ) const { return mom_; }
  G4ThreeVector GetLMom( void ) const { return lmom_; }
  G4int GetTrackNo( void ) const { return trackNo_; }
  G4bool IsTrueSignal() const { return fSignal_; }
  G4double GetXLocal( void ) const { return xl_; }
  G4double GetYLocal( void ) const { return yl_; }
  G4int NumOfHitParticles() const { return HitPartName_.size(); }
  G4String GetHitParticleName( int id ) const 
  { if(id<int(HitPartName_.size())) return HitPartName_[id];
    else return G4String(""); }
  G4double GetPathLength( void ) const { return path_; }
  G4double GetBeta( void ) const { return beta_; }
  G4double GetMass( void ) const { return mass_; }
};

typedef G4THitsCollection<SFTHit> SFTHitsCollection;
extern G4Allocator<SFTHit> SFTHitAllocator;

inline void * SFTHit::operator new( size_t )
{
  return static_cast<void *>( SFTHitAllocator.MallocSingle() );
}

inline void SFTHit::operator delete( void *aHit )
{
  SFTHitAllocator.
    FreeSingle( static_cast<SFTHit *>( aHit ) );
}
		       

#endif
