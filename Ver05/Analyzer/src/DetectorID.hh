/*
 DetectorID.hh
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//DC Number of Plane
const int PlMinBSFT =  1;
const int PlMaxBSFT = 12;

const int PlMinSSFT =  1;
const int PlMaxSSFT = 12;

const int PlMinAFTR = 1;
const int PlMaxAFTR = 2;

const int PlMinAFTL = 1;
const int PlMaxAFTL = 2;

const int PlMinIT1 =  1;
const int PlMaxIT1 = 12;

const int PlMinIT2R = 1;
const int PlMaxIT2R = 9;

const int PlMinIT2L = 1;
const int PlMaxIT2L = 9;

const int PlMinST1 = 1;
const int PlMaxST1 = 9;

const int PlMinST2 = 1;
const int PlMaxST2 = 9;

const int NumOfLayersBSFT = PlMaxBSFT - PlMinBSFT  + 1;
const int NumOfLayersSSFT = PlMaxSSFT - PlMinSSFT  + 1;
const int NumOfLayersAFT = (PlMaxAFTR+PlMaxAFTL) - (PlMinAFTR+PlMinAFTL)  + 2;
const int NumOfLayersIT1 = PlMaxIT1 - PlMinIT1  + 1;
const int NumOfLayersIT2 = (PlMaxIT2R+PlMaxIT2L) - (PlMinIT2R+PlMinIT2L)  + 2;
const int NumOfLayersST1 = PlMaxST1 - PlMinST1  + 1;
const int NumOfLayersST2 = PlMaxST2 - PlMinST2  + 1;

const int NumOfLayersSTIn = (PlMaxSSFT+PlMaxIT1) - (PlMinSSFT+PlMinIT1) + 2;
const int NumOfLayersSTOut = (PlMaxST1+PlMaxST2) - (PlMinST1+PlMinST2) + 2;
const int NumOfLayersITROut = PlMaxIT2R - PlMinIT2R + 1;
const int NumOfLayersITLOut = PlMaxIT2L - PlMinIT2L + 1;

const int NumOfLayersBeamT = NumOfLayersBSFT;
const int NumOfLayersScatInT = NumOfLayersSSFT+NumOfLayersIT1;

const int NumOfLayersScat1T = NumOfLayersSSFT+NumOfLayersIT1+NumOfLayersST1+NumOfLayersST2;
const int NumOfLayersScat2RT = NumOfLayersSSFT+NumOfLayersIT1+NumOfLayersAFT/2+NumOfLayersIT2/2+NumOfLayersST1;
const int NumOfLayersScat2LT = NumOfLayersSSFT+NumOfLayersIT1+NumOfLayersAFT/2+NumOfLayersIT2/2+NumOfLayersST1;



const int PlOffsBSFT = 100;
const int PlOffsSSFT =   0;
const int PlOffsAFT =  30;
const int PlOffsIT1 =  14;
const int PlOffsIT2 =  40;
const int PlOffsST1 =  60;
const int PlOffsST2 =  70;

//Hodo Segments
const int NumOfSegT0  = 1;
const int NumOfSegTof = 1;
const int NumOfSegITofR = 1;
const int NumOfSegITofL = 1;
const int NumOfSegITofSR = 1;
const int NumOfSegITofSL = 1;
const int NumOfSegPADU = 1;
const int NumOfSegPADD = 1;

#endif
