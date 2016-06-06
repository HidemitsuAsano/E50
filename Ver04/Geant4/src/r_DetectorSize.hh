/*
  r_DetectorSize.hh
*/

#ifndef r_DetectorSize_h
#define r_DetectorSize_h 1

/////////////////////
////Real detectors
/////////////////////

//////Tracker size
//BFT
static const G4double rBFT_LaysizeX   = 120.0*mm;
static const G4double rBFT_LaysizeY   = 120.0*mm;
static const G4double rBFT_LaysizeZ   =   1.0*mm;
static const G4double rBFT_TiltAngle  =  30.0*degree;
static const G4double rBFT_BoxThick   =  30.0*mm;
static const G4double rBFT_FrameThick =  10.0*mm;

//SFT
static const G4double rSFT_LaysizeX   =  600.0*mm;
static const G4double rSFT_LaysizeY   =  800.0*mm;
static const G4double rSFT_LaysizeZ   =    0.75*mm;
static const G4double rSFT_TiltAngle  =   30.0*degree;
static const G4double rSFT_FramesizeX = 1800.0*mm;
static const G4double rSFT_BoxThick   =   10.0*mm;
static const G4double rSFT_FrameThick =   30.0*mm;

//AFT
static const G4double rAFT_LaysizeX   =  320.0*mm;
static const G4double rAFT_LaysizeY   = 1200.0*mm;
static const G4double rAFT_LaysizeZ   =    1.0*mm;
static const G4double rAFT_TiltAngle  =    0.0*degree;
static const G4double rAFT_BoxThick   =   10.0*mm;
static const G4double rAFT_FrameThick =   10.0*mm;
static const G4double rAFT_RotAngleR  =  255.0*degree;
static const G4double rAFT_RotAngleL  =  105.0*degree;

//IT1
static const G4double rIT1_LaysizeX     = 1800.0*mm;
static const G4double rIT1_LaysizeUV    = 1800.0*mm;
static const G4double rIT1_LaysizeY     =  800.0*mm;
static const G4double rIT1_LaysizeZ     =    0.01*mm;
static const G4double rIT1_TiltAngle    =    0.0*degree;
static const G4double rIT1_TiltAngle_BT =   45.0*degree;
static const G4double rIT1_BTsize       =   80.0*mm;
static const G4double rIT1_BoxThick     =   30.0*mm;
static const G4double rIT1_FrameThickH  =   50.0*mm;
static const G4double rIT1_FrameThickV  =   80.0*mm;

//IT2
static const G4double rIT2_LaysizeX     = 1300.0*mm;
static const G4double rIT2_LaysizeUV    = 1300.0*mm;
static const G4double rIT2_LaysizeY     = 1400.0*mm;
static const G4double rIT2_LaysizeZ     =    0.01*mm;
static const G4double rIT2_TiltAngle    =    0.0*degree;
static const G4double rIT2_BoxThick     =   30.0*mm;
static const G4double rIT2_FrameThickH  =   50.0*mm;
static const G4double rIT2_FrameThickV  =  100.0*mm;
static const G4double rIT2_RotAngleR    =  255.0*degree;
static const G4double rIT2_RotAngleL    =  105.0*degree;

//ST1
static const G4double rST1_LaysizeX     = 2000.0*mm;
static const G4double rST1_LaysizeUV    = 2000.0*mm;
static const G4double rST1_LaysizeY     = 1300.0*mm;
static const G4double rST1_LaysizeZ     =    0.01*mm;
static const G4double rST1_TiltAngle    =    0.0*degree;
static const G4double rST1_TiltAngle_BT =   45.0*degree;
static const G4double rST1_BoxThick     =   30.0*mm;
static const G4double rST1_FrameThickH  =   50.0*mm;
static const G4double rST1_FrameThickV  =  100.0*mm;
static const G4double rST1_BTsize       =  120.0*mm;
static const G4double rST1_BTpos        =   40.0*mm;

//ST2
static const G4double rST2_LaysizeX     = 2700.0*mm;
static const G4double rST2_LaysizeUV    = 2700.0*mm;
static const G4double rST2_LaysizeY     = 1800.0*mm;
static const G4double rST2_LaysizeZ     =    0.01*mm;
static const G4double rST2_TiltAngle    =    0.0*degree;
static const G4double rST2_TiltAngle_BT =   45.0*degree;
static const G4double rST2_BoxThick     =   30.0*mm;
static const G4double rST2_FrameThickH  =   50.0*mm;
static const G4double rST2_FrameThickV  =  100.0*mm;
static const G4double rST2_BTsize       =  120.0*mm;
static const G4double rST2_BTpos        =   60.0*mm;

//////Counter size
//T0
static const G4double rT0_SegsizeX   = 125.0*mm;
static const G4double rT0_SegsizeY   = 120.0*mm;
static const G4double rT0_SegsizeZ   =   3.0*mm;
static const G4double rT0_SegSpacing = 125.0*mm;
static const G4int    rT0_SegNum     =   1;

//TOF
static const G4double rTof_SegsizeX   = 2800.0*mm;
static const G4double rTof_SegsizeY   = 2000.0*mm;
static const G4double rTof_SegsizeZ   =   30.0*mm;
static const G4double rTof_SegSpacing = 2800.0*mm;
static const G4int    rTof_SegNum     =    1;

//ITof
static const G4double rITof_SegsizeX   = 2000.0*mm;
static const G4double rITof_SegsizeY   = 1500.0*mm;
static const G4double rITof_SegsizeZ   =   30.0*mm;
static const G4double rITof_SegSpacing = 2000.0*mm;
static const G4int    rITof_SegNum     =    1;
static const G4double rITof_RotAngleR  =  255.0*degree;
static const G4double rITof_RotAngleL  =  105.0*degree;

//ITof side wall
static const G4double rITofS_SegsizeX   =  900.0*mm;
static const G4double rITofS_SegsizeY   = 2000.0*mm;
static const G4double rITofS_SegsizeZ   =   30.0*mm;
static const G4double rITofS_SegSpacing =  750.0*mm;
static const G4int    rITofS_SegNum     =    1;
static const G4double rITofS_RotAngleR  =  315.0*degree;
static const G4double rITofS_RotAngleL  =   45.0*degree;

//PAD
static const G4double rPAD_SegsizeX1  =  200.0*mm;
static const G4double rPAD_SegsizeX2  =  100.0*mm;
static const G4double rPAD_SegsizeZ   =   10.0*mm;
static const G4double rPAD_SegsizeY_A = 1750.0*mm;
static const G4double rPAD_SegsizeY_B = 1700.0*mm;
static const G4double rPAD_SegsizeY_C = 1600.0*mm;
static const G4double rPAD_SegsizeY_D = 1400.0*mm;
static const G4double rPAD_SegsizeY_E = 1250.0*mm;
static const G4double rPAD_Offset_A   =    0.0*mm;
static const G4double rPAD_Offset_B   =   25.0*mm;
static const G4double rPAD_SegSpacing =  200.0*mm;
static const G4int    rPAD_SegNum     =   11;

//RICH
static const G4double rRICH_SegsizeX    = 2000.0*mm;
static const G4double rRICH_SegsizeY    = 2000.0*mm;
static const G4double rRICH_SegsizeZ    =   60.0*mm;
static const G4double rRICH_BoxsizeX    = 4500.0*mm;
static const G4double rRICH_BoxsizeY    = 4500.0*mm;
static const G4double rRICH_BoxsizeZ    = 2500.0*mm;
static const G4double rRICH_BoxsizeZ1   =  400.0*mm;
static const G4double rRICH_FrameThick1 =  100.0*mm;
static const G4double rRICH_FrameThick2 =   30.0*mm;
static const G4double rRICH_FrameAngle  =   35.0*degree;

//PID1
static const G4double rPID1_BoxsizeX    = 1100.0*mm;
static const G4double rPID1_BoxsizeY    = 2000.0*mm;
static const G4double rPID1_BoxsizeZ    =  350.0*mm;
static const G4double rPID1_SegsizeZ    =   60.0*mm;
static const G4double rPID1_FrameThick1 =   20.0*mm;

//PID2
static const G4double rPID2_BoxsizeX    = 1500.0*mm;
static const G4double rPID2_BoxsizeY    = 1600.0*mm;
static const G4double rPID2_BoxsizeZ    =  300.0*mm;
static const G4double rPID2_SegsizeZ    =   60.0*mm;
static const G4double rPID2_FrameThick1 =   20.0*mm;
static const G4double rPID2_RotAngleR  =  255.0*degree;
static const G4double rPID2_RotAngleL  =  105.0*degree;

#endif
