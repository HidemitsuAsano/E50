/*
  r_SFT_B.cc

  Segmented SFT
  -> round Fiber + Frame

  Apr. 27th H.Asano
*/

#include "r_SFT_B.hh"
#include "r_DetectorSize.hh"
#include <cmath>


const G4double LzLayer1=-165.0;
const G4double LzLayer2=-135.0;
const G4double LzLayer3=-105.0;
const G4double LzLayer4= -75.0;
const G4double LzLayer5= -45.0;
const G4double LzLayer6= -15.0;
const G4double LzLayer7=  15.0;
const G4double LzLayer8=  45.0;
const G4double LzLayer9=  75.0;
const G4double LzLayer10=105.0;
const G4double LzLayer11=135.0;
const G4double LzLayer12=165.0;

const G4double LzLayer[rSFT_nLayer]={
LzLayer1,
LzLayer2,
LzLayer3,
LzLayer4,
LzLayer5,
LzLayer6,
LzLayer7,
LzLayer8,
LzLayer9,
LzLayer10,
LzLayer11,
LzLayer12};



const G4double MaxDispl= 0.01; /* 10um */

const G4double Deg2Rad = acos(-1.)/180.;
const G4double Rad2Deg = 180./acos(-1.);

#include "globals.hh"
#include "G4Box.hh"
//#include "G4Para.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4ios.hh"

r_SFT_B::r_SFT_B( const G4String & Cname,
		  G4VPhysicalVolume *pMother,
		  const G4RotationMatrix &rotMat,
		  const G4ThreeVector &gPos1, 
		  const G4ThreeVector &gPos2,
		  const G4ThreeVector &gPos3, 
		  const G4ThreeVector &gPos4,
		  const G4ThreeVector &gPos5,
		  const G4ThreeVector &gPos6,
		  const G4ThreeVector &gPos7,
		  const G4ThreeVector &gPos8,
		  const G4ThreeVector &gPos9,
		  const G4ThreeVector &gPos10,
		  const G4ThreeVector &gPos11,
		  const G4ThreeVector &gPos12,
		  const G4ThreeVector &OffsetLocal,
		  G4int id1, G4int id2, G4int id3, 
		  G4int id4, G4int id5, G4int id6,
		  G4int id7, G4int id8, G4int id9,
		  G4int id10, G4int id11, G4int id12,
		  G4Material *matScin, 
		  G4Material *matFrame,
		  G4Material *matArea,
		  G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos6+gPos7)/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double Radius     = rSFT_B_Radius;
  G4double Length     = rSFT_B_Length;

  G4double tilted = rSFT_TiltAngle;
  G4double BaseThick = rSFT_BoxThick*12.0 + 20.0*mm;

  G4double radius_core_in = 0*mm;
  G4double radius_core_out = Radius*rSFT_B_Corefraction;
  G4double radius_innerclad_in = radius_core_out ;
  G4double radius_innerclad_out = radius_innerclad_in+Radius*rSFT_B_InnerCladding;
  G4double radius_outerclad_in = radius_innerclad_out ;
  G4double radius_outerclad_out =  radius_outerclad_in+Radius*rSFT_B_OuterCladding;
  G4double FrameSize = rSFT_FramesizeX;
  G4double FrameT = rSFT_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (rSFT_B_Length+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, rSFT_B_Length/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., rSFT_B_Length/2., BaseThick/2.+0.001*mm );
  

  G4Tubs *solidLayer_core = 
    new G4Tubs( Cname_+"Layer", radius_core_in, radius_core_out, rSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_innerclad = 
    new G4Tubs( Cname_+"Layer_innerclad", radius_innerclad_in, radius_innerclad_out, rSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4Tubs *solidLayer_outerclad = 
    new G4Tubs( Cname_+"Layer_outerclad", radius_outerclad_in, radius_outerclad_out, rSFT_B_Length/2.,0.*deg, 360.*deg);
  
  G4SubtractionSolid *solidFrame =
    new G4SubtractionSolid( Cname_+"Frame", solidBase, solidArea );

   
  logBox   = 
    new G4LogicalVolume( solidBox,  matBox, Cname_+"Box",   0, 0, 0 );
  logFrame   = 
    new G4LogicalVolume( solidFrame,  matFrame, Cname_+"Frame",   0, 0, 0 );
  logArea   = 
    new G4LogicalVolume( solidArea,  matArea, Cname_+"Area",   0, 0, 0 );
  
  logLayer_core = 
    new G4LogicalVolume( solidLayer_core, matScin, Cname_+"Layer_core", 0, 0, 0 );
  
  logLayer_innerclad = 
    new G4LogicalVolume( solidLayer_innerclad, matScin, Cname_+"Layer_innerclad", 0, 0, 0 );
  
  logLayer_outerclad = 
    new G4LogicalVolume( solidLayer_outerclad, matScin, Cname_+"Layer_outerclad", 0, 0, 0 );



  G4VPhysicalVolume *physBox =
    new G4PVPlacement( G4Transform3D( rotMat,
				      G4ThreeVector( DCgPos.x()*mm,
						     DCgPos.y()*mm,
						     DCgPos.z()*mm ) ),
		       Cname_+"Box", logBox, pMother, false, 0 );//name, logical volume, mother volume, no booliean operations its copy number

  G4VPhysicalVolume *physFrame =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Frame", logFrame, physBox, false, 0 );

  G4VPhysicalVolume *physArea =
    new G4PVPlacement( 0, G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm ),
                       Cname_+"Area", logArea, physBox, false, 0 );
 
  // Consitency should be checked.
  G4RotationMatrix InvMat=rotMat.inverse();
  G4ThreeVector lPos1=InvMat*(gPos1-gPosCent);
  G4ThreeVector lPos2=InvMat*(gPos2-gPosCent);
  G4ThreeVector lPos3=InvMat*(gPos3-gPosCent);
  G4ThreeVector lPos4=InvMat*(gPos4-gPosCent);
  G4ThreeVector lPos5=InvMat*(gPos5-gPosCent);
  G4ThreeVector lPos6=InvMat*(gPos6-gPosCent);
  G4ThreeVector lPos7=InvMat*(gPos7-gPosCent);
  G4ThreeVector lPos8=InvMat*(gPos8-gPosCent);
  G4ThreeVector lPos9=InvMat*(gPos9-gPosCent);
  G4ThreeVector lPos10=InvMat*(gPos10-gPosCent);
  G4ThreeVector lPos11=InvMat*(gPos11-gPosCent);
  G4ThreeVector lPos12=InvMat*(gPos12-gPosCent);
  
  if( fabs(lPos1.z()-LzLayer1)>MaxDispl ){
    G4cout << Cname_ << "Layer1: Geometry Error!:" 
	   << "(0,0," << LzLayer1 << ") --> " << lPos1 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos2.z()-LzLayer2)>MaxDispl ){
    G4cout << Cname_ << "Layer2: Geometry Error!:" 
	   << "(0,0," << LzLayer2 << ") --> " << lPos2 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos3.z()-LzLayer3)>MaxDispl ){
    G4cout << Cname_ << "Layer3: Geometry Error!:" 
	   << "(0,0," << LzLayer3 << ") --> " << lPos3 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos4.z()-LzLayer4)>MaxDispl ){
    G4cout << Cname_ << "Layer4: Geometry Error!:" 
	   << "(0,0," << LzLayer4 << ") --> " << lPos4 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos5.z()-LzLayer5)>MaxDispl ){
    G4cout << Cname_ << "Layer5: Geometry Error!:" 
	   << "(0,0," << LzLayer5 << ") --> " << lPos5 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos6.z()-LzLayer6)>MaxDispl ){
    G4cout << Cname_ << "Layer6: Geometry Error!:" 
	   << "(0,0," << LzLayer6 << ") --> " << lPos6 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos7.z()-LzLayer7)>MaxDispl ){
    G4cout << Cname_ << "Layer7: Geometry Error!:" 
	   << "(0,0," << LzLayer7 << ") --> " << lPos7 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos8.z()-LzLayer8)>MaxDispl ){
    G4cout << Cname_ << "Layer8: Geometry Error!:" 
	   << "(0,0," << LzLayer8 << ") --> " << lPos8 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos9.z()-LzLayer9)>MaxDispl ){
    G4cout << Cname_ << "Layer9: Geometry Error!:" 
	   << "(0,0," << LzLayer9 << ") --> " << lPos9 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos10.z()-LzLayer10)>MaxDispl ){
    G4cout << Cname_ << "Layer10: Geometry Error!:" 
	   << "(0,0," << LzLayer10 << ") --> " << lPos10 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos11.z()-LzLayer11)>MaxDispl ){
    G4cout << Cname_ << "Layer11: Geometry Error!:" 
	   << "(0,0," << LzLayer11 << ") --> " << lPos11 << G4endl;
    G4Exception("Geometry Error");
  }
  if( fabs(lPos12.z()-LzLayer12)>MaxDispl ){
    G4cout << Cname_ << "Layer12: Geometry Error!:" 
	   << "(0,0," << LzLayer12 << ") --> " << lPos12 << G4endl;
    G4Exception("Geometry Error");
  }

  const G4double offset_2ndlayer = 0.5*mm;
  for( int id=1; id<=rSFT_SegNum; ++id ){ 
    G4double ofsScinX=-rSFT_SegSpacing*(rSFT_SegNum/2.+0.5-id);
     
    
    // Layer1-1 x
    G4RotationMatrix *layrot = new G4RotationMatrix();
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer1*mm ),
    		       Cname_+"Layer1_1_core", logLayer_core, physArea, false, 1000*id1+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer1*mm ),
    		       Cname_+"Layer1_1_innerclad", logLayer_innerclad, physArea, false, 1000*id1+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer1*mm ),
    		       Cname_+"Layer1_1_outerclad", logLayer_outerclad, physArea, false, 1000*id1+id );

    
    // Layer1-2 x
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_core", logLayer_core, physArea, false, 1000*id1+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_innerclad", logLayer_innerclad, physArea, false, 1000*id1+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_outerclad", logLayer_outerclad, physArea, false, 1000*id1+id+n_SFT_B_fiber);
    
    
    // Layer2-1 u
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (-1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer2*mm ),
    		       Cname_+"Layer2_1_core", logLayer_core, physArea, false, 1000*id2+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer2*mm ),
    		       Cname_+"Layer2_1_innerclad", logLayer_innerclad, physArea, false, 1000*id2+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer2*mm ),
    		       Cname_+"Layer2_1_outerclad", logLayer_outerclad, physArea, false, 1000*id2+id );

    // Layer2-2 u
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer2*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer2_2_core", logLayer_core, physArea, false, 1000*id2+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer2*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer2_2_innerclad", logLayer_innerclad, physArea, false, 1000*id2+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer2*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer2_2_outerclad", logLayer_outerclad, physArea, false, 1000*id2+id+n_SFT_B_fiber);
    
    
    // Layer3-1 v
    
    layrot = new G4RotationMatrix();
   
    layrot->rotateZ( (1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer3*mm ),
    		       Cname_+"Layer3_1_core", logLayer_core, physArea, false, 1000*id3+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer3*mm ),
    		       Cname_+"Layer3_1_innerclad", logLayer_innerclad, physArea, false, 1000*id3+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer3*mm ),
    		       Cname_+"Layer3_1_outerclad", logLayer_outerclad, physArea, false, 1000*id3+id );

    // Layer3-2 v
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer3*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer3_2_core", logLayer_core, physArea, false, 1000*id3+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer3*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer3_2_innerclad", logLayer_innerclad, physArea, false, 1000*id3+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer3*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer3_2_outerclad", logLayer_outerclad, physArea, false, 1000*id3+id+n_SFT_B_fiber);
    
    

    
    // Layer4-1 x
    layrot = new G4RotationMatrix();
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer4*mm ),
    		       Cname_+"Layer4_1_core", logLayer_core, physArea, false, 1000*id4+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer4*mm ),
    		       Cname_+"Layer4_1_innerclad", logLayer_innerclad, physArea, false, 1000*id4+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer4*mm ),
    		       Cname_+"Layer4_1_outerclad", logLayer_outerclad, physArea, false, 1000*id4+id );

    // Layer4-2 x
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer4*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer4_2_core", logLayer_core, physArea, false, 1000*id4+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer4*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer4_2_innerclad", logLayer_innerclad, physArea, false, 1000*id4+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer4*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer4_2_outerclad", logLayer_outerclad, physArea, false, 1000*id4+id+n_SFT_B_fiber);
    
    
    // Layer5-1 u
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (-1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer5*mm ),
    		       Cname_+"Layer5_1_core", logLayer_core, physArea, false, 1000*id5+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer5*mm ),
    		       Cname_+"Layer5_1_innerclad", logLayer_innerclad, physArea, false, 1000*id5+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer5*mm ),
    		       Cname_+"Layer5_1_outerclad", logLayer_outerclad, physArea, false, 1000*id5+id );

    // Layer5-2 u
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer5*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer5_2_core", logLayer_core, physArea, false, 1000*id5+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer5*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer5_2_innerclad", logLayer_innerclad, physArea, false, 1000*id5+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer5*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer5_2_outerclad", logLayer_outerclad, physArea, false, 1000*id5+id+n_SFT_B_fiber);
    
    
    
    // Layer6-1 v
    
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer6*mm ),
    		       Cname_+"Layer6_1_core", logLayer_core, physArea, false, 1000*id6+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer6*mm ),
    		       Cname_+"Layer6_1_innerclad", logLayer_innerclad, physArea, false, 1000*id6+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer6*mm ),
    		       Cname_+"Layer6_1_outerclad", logLayer_outerclad, physArea, false, 1000*id6+id );

    // Layer6-2
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer6*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer6_2_core", logLayer_core, physArea, false, 1000*id6+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer6*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer6_2_innerclad", logLayer_innerclad, physArea, false, 1000*id6+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer6*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer6_2_outerclad", logLayer_outerclad, physArea, false, 1000*id6+id+n_SFT_B_fiber);
    
    
    // Layer7-1 u
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (-1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer7*mm ),
    		       Cname_+"Layer7_1_core", logLayer_core, physArea, false, 1000*id7+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer7*mm ),
    		       Cname_+"Layer7_1_innerclad", logLayer_innerclad, physArea, false, 1000*id7+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer7*mm ),
    		       Cname_+"Layer7_1_outerclad", logLayer_outerclad, physArea, false, 1000*id7+id );

    
    // Layer7-2 u
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer7*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer7_2_core", logLayer_core, physArea, false, 1000*id7+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer7*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer7_2_innerclad", logLayer_innerclad, physArea, false, 1000*id7+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer7*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer7_2_outerclad", logLayer_outerclad, physArea, false, 1000*id7+id+n_SFT_B_fiber);
    
    
    // Layer8-1 v
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer8*mm ),
    		       Cname_+"Layer8_1_core", logLayer_core, physArea, false, 1000*id8+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer8*mm ),
    		       Cname_+"Layer8_1_innerclad", logLayer_innerclad, physArea, false, 1000*id8+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer8*mm ),
    		       Cname_+"Layer8_1_outerclad", logLayer_outerclad, physArea, false, 1000*id8+id );

    // Layer8-2 v
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer8*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer8_2_core", logLayer_core, physArea, false, 1000*id8+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer8*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer8_2_innerclad", logLayer_innerclad, physArea, false, 1000*id8+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer8*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer8_2_outerclad", logLayer_outerclad, physArea, false, 1000*id8+id+n_SFT_B_fiber);
    
    
    // Layer9-1 x
    layrot = new G4RotationMatrix();
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer9*mm ),
    		       Cname_+"Layer9_1_core", logLayer_core, physArea, false, 1000*id9+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer9*mm ),
    		       Cname_+"Layer9_1_innerclad", logLayer_innerclad, physArea, false, 1000*id9+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer9*mm ),
    		       Cname_+"Layer9_1_outerclad", logLayer_outerclad, physArea, false, 1000*id9+id );

    // Layer9-2 x
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer9*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer9_2_core", logLayer_core, physArea, false, 1000*id9+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer9*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer9_2_innerclad", logLayer_innerclad, physArea, false, 1000*id9+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer9*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer9_2_outerclad", logLayer_outerclad, physArea, false, 1000*id9+id+n_SFT_B_fiber);
    
    

    // Layer10-1 u
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (-1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer10*mm ),
    		       Cname_+"Layer10_1_core", logLayer_core, physArea, false, 1000*id10+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer10*mm ),
    		       Cname_+"Layer10_1_innerclad", logLayer_innerclad, physArea, false, 1000*id10+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer10*mm ),
    		       Cname_+"Layer10_1_outerclad", logLayer_outerclad, physArea, false, 1000*id10+id );

    // Layer10-2 u
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer10*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer10_2_core", logLayer_core, physArea, false, 1000*id10+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer10*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer10_2_innerclad", logLayer_innerclad, physArea, false, 1000*id10+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer10*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer10_2_outerclad", logLayer_outerclad, physArea, false, 1000*id10+id+n_SFT_B_fiber);
    
    
    // Layer11-1 v
    layrot = new G4RotationMatrix();
    layrot->rotateZ( (1.)*rSFT_TiltAngle);
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer11*mm ),
    		       Cname_+"Layer11_1_core", logLayer_core, physArea, false, 1000*id11+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer11*mm ),
    		       Cname_+"Layer11_1_innerclad", logLayer_innerclad, physArea, false, 1000*id11+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer11*mm ),
    		       Cname_+"Layer11_1_outerclad", logLayer_outerclad, physArea, false, 1000*id11+id );

    // Layer11-2 v
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer11*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer11_2_core", logLayer_core, physArea, false, 1000*id11+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer11*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer11_2_innerclad", logLayer_innerclad, physArea, false, 1000*id11+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer11*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer11_2_outerclad", logLayer_outerclad, physArea, false, 1000*id11+id+n_SFT_B_fiber);
    
    
    
    // Layer12-1 x
    layrot = new G4RotationMatrix();
    layrot->rotateX( 90.*deg);
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer12*mm ),
    		       Cname_+"Layer12_1_core", logLayer_core, physArea, false, 1000*id12+id );

    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer12*mm ),
    		       Cname_+"Layer12_1_innerclad", logLayer_innerclad, physArea, false, 1000*id12+id );
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer12*mm ),
    		       Cname_+"Layer12_1_outerclad", logLayer_outerclad, physArea, false, 1000*id12+id );

    // Layer12-2
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer12*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer12_2_core", logLayer_core, physArea, false, 1000*id12+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer12*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer12_2_innerclad", logLayer_innerclad, physArea, false, 1000*id12+id+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer12*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer12_2_outerclad", logLayer_outerclad, physArea, false, 1000*id12+id+n_SFT_B_fiber);
    
    
    
  }
}

void r_SFT_B::SetVisAttributes( const G4VisAttributes *attLayer,
        const G4VisAttributes *attInnerClad,
        const G4VisAttributes *attOuterClad,
				const G4VisAttributes *attFrame,
				const G4VisAttributes *attArea,
				const G4VisAttributes *attBox )
{
  if(attLayer) logLayer_core->SetVisAttributes(attLayer); 
  if(attLayer) logLayer_innerclad->SetVisAttributes(attInnerClad); 
  if(attLayer) logLayer_outerclad->SetVisAttributes(attOuterClad); 
  if(attFrame) logFrame->SetVisAttributes(attFrame); 
  if(attArea)  logArea->SetVisAttributes(attArea); 
  if(attBox)   logBox->SetVisAttributes(attBox); 
}

void r_SFT_B::SetSensitiveDetector( G4VSensitiveDetector *SD )
{
  logLayer_core->SetSensitiveDetector(SD);
}
