/*
  r_DetectorSize.hh
*/

#ifndef r_DetectorSize_h
#define r_DetectorSize_h 1

/////////////////////
////Simple detectors
/////////////////////

//////Tracker size
//SFT
static const G4double rSFT_LaysizeX   =    1.0*mm;
static const G4double rSFT_LaysizeY   =  200.0*mm;
static const G4double rSFT_LaysizeZ   =    1.0*mm;
static const G4double rSFT_TiltAngle  =   30.0*degree;
static const G4double rSFT_SegSpacing =    1.0*mm;
static const G4int    rSFT_SegNum     =   10;
static const G4double rSFT_FramesizeX =  200.0*mm;
static const G4double rSFT_BoxThick   =   30.0*mm;
static const G4double rSFT_FrameThick =   30.0*mm;

static const G4double rSFT_livefraction = 0.96;//fraction of senstive area of fiber

//////Counter size
//T0
static const G4double rT0_SegsizeX   =   9.0*mm;
static const G4double rT0_SegsizeY   = 100.0*mm;
static const G4double rT0_SegsizeZ   =   3.0*mm;
static const G4double rT0_SegSpacing =   3.1*mm;
static const G4int    rT0_SegNum     =   1;

#endif
