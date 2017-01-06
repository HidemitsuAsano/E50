/*
  r_DetectorSize.hh
*/

#ifndef r_DetectorSize_h
#define r_DetectorSize_h 1

/////////////////////
////Realistic detectors
/////////////////////

static const G4int    rSFT_nLayer = 12;
//////Tracker size
//SFT 1mm square fiber Single Cladding
static const G4double rSFT_Square_LaysizeX   =    1.0*mm;
static const G4double rSFT_Square_LaysizeY   =  200.0*mm;
static const G4double rSFT_Square_LaysizeZ   =    1.0*mm;
static const G4double rSFT_SegSpacing =    1.0*mm;
static const G4int    rSFT_SegNum     =   100;
static const G4double rSFT_FramesizeX =  200.0*mm;
static const G4double rSFT_BoxThick   =   30.0*mm;
static const G4double rSFT_FrameThick =   30.0*mm;

static const G4double rSFT_livefraction = 0.96;//fraction of senstive area of fiber

//SFT round fiber Multi Cladding phi=1mm
static const G4double rSFT_Round_Radius     =    0.5*mm;
static const G4double rSFT_Round_Length     =  200.0*mm;
static const G4double rSFT_Round_SegSpacing =    1.00*mm;//must be same as dWdX in param file
//static const G4int    rSFT_Round_SegNum     =   10;//# prototype
static const G4int    rSFT_Round_SegNum     =   100;//# of fiber (or bundle) per one sublayer
static const G4double rSFT_Round_FramesizeX =  200.0*mm;
static const G4double rSFT_Round_BoxThick   =   30.0*mm;
static const G4double rSFT_Round_FrameThick =   30.0*mm;

static const G4double rSFT_Round_Corefraction = 0.94;//fraction of senstive area of fiber
static const G4double rSFT_Round_InnerCladding = 0.03;//fraction of inner cladding of fiber
static const G4double rSFT_Round_OuterCladding = 0.03;//fraction of outer cladding of fiber

//SFT_C 0.5mm square fiber Single Cladding
//static const G4double rSFT_C_LaysizeX   =    0.2*mm;
//static const G4double rSFT_C_LaysizeY   =  200.0*mm;
//static const G4double rSFT_C_LaysizeZ   =    0.2*mm;
//static const G4double rSFT_C_SegSpacing =    0.2*mm;
//static const G4int    rSFT_C_SegNum     =   50;


//SFT_D 0.5 mm round fibe Multi Cladding r=0.25*mm
static const G4double rSFT_D_Radius     =    0.25*mm;
static const G4double rSFT_D_Length     =  200.0*mm;
static const G4double rSFT_D_TiltAngle  =   30.0*degree;
static const G4double rSFT_D_SegSpacing =    1.0*mm;
static const G4int    rSFT_D_SegNum     =   20;
static const G4double rSFT_D_FramesizeX =  200.0*mm;
static const G4double rSFT_D_BoxThick   =   30.0*mm;
static const G4double rSFT_D_FrameThick =   30.0*mm;

static const G4double rSFT_D_Corefraction = 0.94;//fraction of senstive area of fiber
static const G4double rSFT_D_InnerCladding = 0.03;//fraction of inner cladding of fiber
static const G4double rSFT_D_OuterCladding = 0.03;//fraction of outer cladding of fiber
//////Counter size
//T0
static const G4double rT0_SegsizeX   =   9.0*mm;
static const G4double rT0_SegsizeY   = 100.0*mm;
static const G4double rT0_SegsizeZ   =   3.0*mm;
static const G4double rT0_SegSpacing =   3.1*mm;
static const G4int    rT0_SegNum     =   1;

#endif
