/*
  r_SFT_B.cc

  Segmented SFT
  -> round Fiber + Frame

  Apr. 27th H.Asano
*/

#include "r_SFT_B.hh"
#include "r_DetectorSize.hh"
#include <cmath>
#include <sstream>
#include <string>

const G4double LzLayer[rSFT_nLayer]={
-165.0*mm,
-135.0*mm,
-105.0*mm,
 -75.0*mm,
 -45.0*mm,
 -15.0*mm,
  15.0*mm,
  45.0*mm,
  75.0*mm,
 105.0*mm,
 135.0*mm,
 165.0*mm
};



const G4double MaxDispl= 0.01; /* 10um */


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
		  const G4ThreeVector *gPos, 
		  const G4ThreeVector &OffsetLocal,
		  G4int *detid,
		  G4Material *matScin, 
		  G4Material *matPMMA, 
		  G4Material *matFP, 
		  G4Material *matFrame,
		  G4Material *matArea,
		  G4Material *matBox )
  : Cname_(Cname)
{
  G4ThreeVector gPosCent=(gPos[5]+gPos[6])/2.;
  G4ThreeVector gOffset=rotMat*OffsetLocal;
  G4ThreeVector DCgPos=gPosCent+gOffset;
  
  G4double Radius     = rSFT_B_Radius;

  G4double BaseThick = rSFT_BoxThick*12.0 + 20.0*mm;

  G4double radius_core_in = 0*mm;
  G4double radius_core_out = Radius*rSFT_B_Corefraction;
  G4double radius_innerclad_in = radius_core_out ;
  G4double radius_innerclad_out = radius_innerclad_in+Radius*rSFT_B_InnerCladding;
  G4double radius_outerclad_in = radius_innerclad_out ;
  G4double radius_outerclad_out =  radius_outerclad_in+Radius*rSFT_B_OuterCladding;
  G4double FrameSize = rSFT_B_FramesizeX;
  G4double FrameT = rSFT_B_FrameThick;

  G4Box *solidBox  = 
    new G4Box( Cname_+"Box", (FrameSize+1.0*mm)/2.+FrameT, (rSFT_B_Length+1.0*mm)/2.+FrameT, (BaseThick+0.5*mm)/2. );
  G4Box *solidBase  = 
    new G4Box( Cname_+"Base", FrameSize/2.+FrameT, rSFT_B_Length/2.+FrameT, BaseThick/2. );
  G4Box *solidArea  = 
    new G4Box( Cname_+"Area", FrameSize/2., rSFT_B_Length/2., BaseThick/2.+0.001*mm );
  //   new G4Box( Cname_+"Area", FrameSize/2., rSFT_B_Length/2., BaseThick/2. );
  

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
    new G4LogicalVolume( solidLayer_innerclad, matPMMA, Cname_+"Layer_innerclad", 0, 0, 0 );
  
  logLayer_outerclad = 
    new G4LogicalVolume( solidLayer_outerclad, matFP, Cname_+"Layer_outerclad", 0, 0, 0 );



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
  G4ThreeVector lPos[12];
  for(G4int ilr = 0; ilr < rSFT_nLayer;ilr++){
    lPos[ilr] = InvMat*(gPos[ilr]-gPosCent);
  
    if( fabs(lPos[ilr].z()-LzLayer[ilr])>MaxDispl){
      G4cout << Cname_ << "Layer" << ilr+1 << ": Geometry Error!:" 
        << "(0,0," << LzLayer[ilr] << ") --> " << lPos[ilr] << G4endl;
      G4Exception("Geometry Error");
    }
  }
 

  const G4double offset_2ndsublayer = 0.5*mm;
  //x u v x u v u v x u v x 
  //Here
  //x: 0 
  //u: 1
  //v: 2
  const G4int layer_conf[rSFT_nLayer]={0, 1, 2, 0, 1,2 ,1,2,0,1,2,0}; 
  
  G4RotationMatrix *layrot[3];
  //x layer local -> global rotation
  layrot[0] = new G4RotationMatrix();
  layrot[0]->rotateX(90.*deg);
  
  //u layer local ->global rotation
  layrot[1] = new G4RotationMatrix();
  layrot[1]->rotateZ( (-1.000)*rSFT_TiltAngle);
  layrot[1]->rotateX( 90.*deg);
  
  //v layer local ->global rotation
  layrot[2] = new G4RotationMatrix();
  layrot[2]->rotateZ( (1.000)*rSFT_TiltAngle);
  layrot[2]->rotateX( 90.*deg);

    
    for( G4int iseg=1; iseg<=rSFT_SegNum; ++iseg ){ 
      G4double ofsScinX=-rSFT_SegSpacing*(rSFT_SegNum/2.+0.5-iseg);
      for( G4int ilr=0; ilr<rSFT_nLayer;ilr++){
        G4int conf = layer_conf[ilr];
      
        //set name
        std::ostringstream os;
        os << ilr+1;
        G4String Lname = os.str();;//+"_"+1;

        //first sublayer
        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_core", logLayer_core, physArea, false, 1000*detid[ilr]+iseg );

        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+iseg );

        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX, 0.0*mm, -LzLayer[ilr] ),
            Cname_+"Layer"+Lname+"_1_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+iseg );

        //second sublayer
        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*1.732051 ),
  //      new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*1.73205 ),
     //   new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*2. ),
            Cname_+"Layer"+Lname+"_2_core", logLayer_core, physArea, false, 1000*detid[ilr]+iseg+n_SFT_B_fiber );

        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*1.732051 ),
     //   new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*2. ),
            Cname_+"Layer"+Lname+"_2_innerclad", logLayer_innerclad, physArea, false, 1000*detid[ilr]+iseg+n_SFT_B_fiber );

        new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*1.732051 ),
     //   new G4PVPlacement( layrot[conf], G4ThreeVector(  ofsScinX+offset_2ndsublayer, 0.0*mm, -LzLayer[ilr]+rSFT_B_Radius*2.0 ),
            Cname_+"Layer"+Lname+"_2_outerclad", logLayer_outerclad, physArea, false, 1000*detid[ilr]+iseg+n_SFT_B_fiber );
      }
    /*
    // Layer1-2 x
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_core", logLayer_core, physArea, false, 1000*id1+iseg+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_innerclad", logLayer_innerclad, physArea, false, 1000*id1+iseg+n_SFT_B_fiber);
    
    new G4PVPlacement( layrot, G4ThreeVector(  ofsScinX+offset_2ndlayer, 0.0*mm, -LzLayer1*mm+rSFT_B_Radius*1.732051 ),
    		       Cname_+"Layer1_2_outerclad", logLayer_outerclad, physArea, false, 1000*id1+iseg+n_SFT_B_fiber);
    
    
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
    
    */    
    
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
