/*
  s_DetectorSize.hh
*/

#ifndef s_DetectorSize_h
#define s_DetectorSize_h 1

/////////////////////
////Simple detectors
/////////////////////

//////Tracker size
//BFT
static const G4double BFT_LaysizeX   = 120.0*mm;
static const G4double BFT_LaysizeY   = 120.0*mm;
static const G4double BFT_LaysizeZ   =   1.0*mm;
static const G4double BFT_TiltAngle  =  30.0*degree;
static const G4double BFT_BoxThick   =  30.0*mm;
static const G4double BFT_FrameThick =  10.0*mm;

//SFT
static const G4double SFT_LaysizeX   =  600.0*mm;
static const G4double SFT_LaysizeY   =  800.0*mm;
static const G4double SFT_LaysizeZ   =    0.75*mm;
static const G4double SFT_TiltAngle  =   30.0*degree;
static const G4double SFT_FramesizeX = 1800.0*mm;
static const G4double SFT_BoxThick   =   10.0*mm;
static const G4double SFT_FrameThick =   30.0*mm;

//AFT
static const G4double AFT_LaysizeX   =  320.0*mm;
static const G4double AFT_LaysizeY   = 1200.0*mm;
static const G4double AFT_LaysizeZ   =    1.0*mm;
static const G4double AFT_TiltAngle  =    0.0*degree;
static const G4double AFT_BoxThick   =   10.0*mm;
static const G4double AFT_FrameThick =   10.0*mm;

//IT1
static const G4double IT1_LaysizeX     = 1800.0*mm;
static const G4double IT1_LaysizeUV    = 1800.0*mm;
static const G4double IT1_LaysizeY     =  800.0*mm;
static const G4double IT1_LaysizeZ     =    0.01*mm;
static const G4double IT1_TiltAngle    =    0.0*degree;
static const G4double IT1_TiltAngle_BT =   45.0*degree;
static const G4double IT1_BTsize       =   80.0*mm;
static const G4double IT1_BoxThick     =   30.0*mm;
static const G4double IT1_FrameThickH  =   50.0*mm;
static const G4double IT1_FrameThickV  =   80.0*mm;

//IT2
static const G4double IT2_LaysizeX     = 1500.0*mm;
static const G4double IT2_LaysizeUV    = 1500.0*mm;
static const G4double IT2_LaysizeY     = 1400.0*mm;
static const G4double IT2_LaysizeZ     =    0.01*mm;
static const G4double IT2_TiltAngle    =    0.0*degree;
static const G4double IT2_BoxThick     =   30.0*mm;
static const G4double IT2_FrameThickH  =   50.0*mm;
static const G4double IT2_FrameThickV  =  100.0*mm;

//ST1
static const G4double ST1_LaysizeX     = 2000.0*mm;
static const G4double ST1_LaysizeUV    = 2000.0*mm;
static const G4double ST1_LaysizeY     = 1300.0*mm;
static const G4double ST1_LaysizeZ     =    0.01*mm;
static const G4double ST1_TiltAngle    =    0.0*degree;
static const G4double ST1_TiltAngle_BT =   45.0*degree;
static const G4double ST1_BoxThick     =   30.0*mm;
static const G4double ST1_FrameThickH  =   50.0*mm;
static const G4double ST1_FrameThickV  =  100.0*mm;
static const G4double ST1_BTsize       =  120.0*mm;
static const G4double ST1_BTpos        =   40.0*mm;

//ST2
static const G4double ST2_LaysizeX     = 2700.0*mm;
static const G4double ST2_LaysizeUV    = 2700.0*mm;
static const G4double ST2_LaysizeY     = 1800.0*mm;
static const G4double ST2_LaysizeZ     =    0.01*mm;
static const G4double ST2_TiltAngle    =    0.0*degree;
static const G4double ST2_TiltAngle_BT =   45.0*degree;
static const G4double ST2_BoxThick     =   30.0*mm;
static const G4double ST2_FrameThickH  =   50.0*mm;
static const G4double ST2_FrameThickV  =  100.0*mm;
static const G4double ST2_BTsize       =  120.0*mm;
static const G4double ST2_BTpos        =   60.0*mm;

//////Counter size
//T0
static const G4double T0_SegsizeX   = 125.0*mm;
static const G4double T0_SegsizeY   = 120.0*mm;
static const G4double T0_SegsizeZ   =   3.0*mm;
static const G4double T0_SegSpacing = 125.0*mm;
static const G4int    T0_SegNum     =   1;

//TOF
static const G4double Tof_SegsizeX   = 2800.0*mm;
static const G4double Tof_SegsizeY   = 2000.0*mm;
static const G4double Tof_SegsizeZ   =   30.0*mm;
static const G4double Tof_SegSpacing = 2800.0*mm;
static const G4int    Tof_SegNum     =    1;

//ITof
static const G4double ITof_SegsizeX   = 2300.0*mm;
static const G4double ITof_SegsizeY   = 1500.0*mm;
static const G4double ITof_SegsizeZ   =   30.0*mm;
static const G4double ITof_SegSpacing = 2000.0*mm;
static const G4int    ITof_SegNum     =    1;
static const G4double ITof_RotAngleR  =  255.0*degree;
static const G4double ITof_RotAngleL  =  105.0*degree;

//ITof side wall
static const G4double ITofS_SegsizeX   =  700.0*mm;
static const G4double ITofS_SegsizeY   = 2000.0*mm;
static const G4double ITofS_SegsizeZ   =   30.0*mm;
static const G4double ITofS_SegSpacing =  750.0*mm;
static const G4int    ITofS_SegNum     =    1;
static const G4double ITofS_RotAngleR  =  330.0*degree;
static const G4double ITofS_RotAngleL  =   30.0*degree;

//PAD
static const G4double PAD_SegsizeX1  =  200.0*mm;
static const G4double PAD_SegsizeX2  =  100.0*mm;
static const G4double PAD_SegsizeZ   =   10.0*mm;
static const G4double PAD_SegsizeY_A = 1750.0*mm;
static const G4double PAD_SegsizeY_B = 1700.0*mm;
static const G4double PAD_SegsizeY_C = 1600.0*mm;
static const G4double PAD_SegsizeY_D = 1400.0*mm;
static const G4double PAD_SegsizeY_E = 1250.0*mm;
static const G4double PAD_Offset_A   =    0.0*mm;
static const G4double PAD_Offset_B   =   25.0*mm;
static const G4double PAD_SegSpacing =  200.0*mm;
static const G4int    PAD_SegNum     =   11;

//RICH
static const G4double RICH_SegsizeX    = 2000.0*mm;
static const G4double RICH_SegsizeY    = 2000.0*mm;
static const G4double RICH_SegsizeZ    =   60.0*mm;
static const G4double RICH_BoxsizeX    = 4500.0*mm;
static const G4double RICH_BoxsizeY    = 4500.0*mm;
static const G4double RICH_BoxsizeZ    = 2500.0*mm;
static const G4double RICH_BoxsizeZ1   =  400.0*mm;
static const G4double RICH_FrameThick1 =  100.0*mm;
static const G4double RICH_FrameThick2 =   30.0*mm;
static const G4double RICH_FrameAngle  =   35.0*degree;

//PID1
static const G4double PID1_BoxsizeX    = 1100.0*mm;
static const G4double PID1_BoxsizeY    = 2000.0*mm;
static const G4double PID1_BoxsizeZ    =  350.0*mm;
static const G4double PID1_SegsizeZ    =   60.0*mm;
static const G4double PID1_FrameThick1 =   20.0*mm;

//PID2
static const G4double PID2_BoxsizeX    = 1700.0*mm;
static const G4double PID2_BoxsizeY    = 1600.0*mm;
static const G4double PID2_BoxsizeZ    =  300.0*mm;
static const G4double PID2_SegsizeZ    =   60.0*mm;
static const G4double PID2_FrameThick1 =   20.0*mm;
static const G4double PID2_RotAngleR  =  255.0*degree;
static const G4double PID2_RotAngleL  =  105.0*degree;

#endif
