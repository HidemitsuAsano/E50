/*
  s_DetectorSize.hh
*/

#ifndef s_DetectorSize_h
#define s_DetectorSize_h 1

/////////////////////
////Simple detectors
/////////////////////

//////Tracker size
//SFT
static const G4double sSFT_LaysizeX   =   10.0*mm;
static const G4double sSFT_LaysizeY   =  200.0*mm;
static const G4double sSFT_LaysizeZ   =    1.0*mm;
static const G4double sSFT_TiltAngle  =   30.0*degree;
static const G4double sSFT_FramesizeX =  200.0*mm;
static const G4double sSFT_BoxThick   =   30.0*mm;
static const G4double sSFT_FrameThick =   30.0*mm;

//////Counter size
//T0
static const G4double sT0_SegsizeX   =   9.0*mm;
static const G4double sT0_SegsizeY   = 100.0*mm;
static const G4double sT0_SegsizeZ   =   3.0*mm;
static const G4double sT0_SegSpacing =   3.1*mm;
static const G4int    sT0_SegNum     =   1;

#endif
