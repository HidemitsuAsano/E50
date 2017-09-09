#include "TrLocalTrack.hh"
//#include "TrLTrackHit.hh"
#include "SFTCluster.hh"
#include "GeomMan.hh"
#include "BSFTReco.hh"
#include "DetectorID.hh"
#include "MathTools.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

const double Deg2Rad = acos(-1.)/180.;
const double Rad2Deg = 180./acos(-1.);

const int ReservedNumOfHits = 16;
//const unsigned int TrLocalMinNHits  = 6;
const unsigned int TrLocalMinNHits  = 4;
const unsigned int TrLocalMinNHitsVXU = 3;
const unsigned int TrLocalMinNHitsVXU2= 4;
const unsigned int TrLocalMinNHits2 = 8;


TrLocalTrack::TrLocalTrack()
  : status_(false), 
    x0_(0.0), y0_(0.0), 
    u0_(0.0), v0_(0.0),
    u1_(0.0), v1_(0.0),
    a_(0.0),b_(0.0),c_(0.0),
    chisqr_(0),
    gftstatus_(true),
    trackid_(-1)
{
  sftclusterArray_.reserve( ReservedNumOfHits );
  Coefficients_x[0]=0;
  Coefficients_y[0]=0;
}

TrLocalTrack::~TrLocalTrack()
{
}

bool TrLocalTrack::allocateBufferArea( void )
{
  double *buf = new double [6*12];
  if(buf==0){
    return false;
  }
  for(int i=0; i<6; ++i ) {
    Coefficients_x[i]= &(buf[i*12]);
    Coefficients_y[i]= &(buf[i*12]);
  }
  return true;
}

SFTCluster *TrLocalTrack::GetHit( std::size_t nth ) const
{
  if( nth<sftclusterArray_.size() )
    return sftclusterArray_[nth];
  else
    return 0;
}


SFTCluster *TrLocalTrack::GetHitOfLayerNumber( int lnum ) const
{
  for( std::size_t i=0; i<sftclusterArray_.size(); ++i )
    if( sftclusterArray_[i]->GetLayer()==lnum )
      return sftclusterArray_[i];
  return 0;
}


//liner fitting
bool TrLocalTrack::DoFit( void )
{
  const std::string funcname = "[TrLocalTrack::DoFit()]";
  //const GeomMan & geomMan=GeomMan::GetInstance();
  
  std::size_t n = sftclusterArray_.size();
  static int state = 0;
  if(!state){
    std::cout << __FILE__ << " L" << __LINE__ << " " <<  __FUNCTION__ <<" Minimum Numbef of hits: "  << TrLocalMinNHits << std::endl;
    state++;
  }

  if(n < TrLocalMinNHits ) return status_ = false;
  
  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);
  
  //bool hitpattern[12]={false};
  for( std::size_t i=0; i<n; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      //double ww = geomMan.GetResolution( lnum );

      //assuming uniform distribution, ADC info is not used here
      //perhaps, fiber diameter/sqrt(12) because left-right ambiguity is reduced by 2 layer structure
      double ww = sftclusterp->GetClusterLxSize()/sqrt(12.0);
      
      double gz = sftclusterp->GetGlobalZ();
      
      //std::cout << __FILE__ << " line " << __LINE__ << 
      //" resolution " << ww << " gz  " << gz << std::endl;
      
      //flag for removing some layers 
      bool isChangeTrackingConf = false;
      static bool isStated = false; 

      int lnum = sftclusterp->GetLayer();
      //hitpattern[lnum] = true;
      if(isChangeTrackingConf){
        //if(lnum == 0) ww = 999999.9;//x
        if(lnum == 1) ww = 999999.9;//u
        if(lnum == 2) ww = 999999.9;//v
        //if(lnum == 3) ww = 999999.9;//x
        if(lnum == 4) ww = 999999.9;//u
        if(lnum == 5) ww = 999999.9;//v
        if(lnum == 6) ww = 999999.9;//u
        if(lnum == 7) ww = 999999.9;//v
        //if(lnum == 8) ww = 999999.9;//x
        if(lnum == 9) ww = 999999.9;//u
        if(lnum == 10) ww = 999999.9;//v
        //if(lnum == 11) ww = 999999.9;//x
        if(!isStated){
          std::cout << __FILE__ << " l. " << __LINE__ << "some layers are removed from track fitting" << std::endl;
          isStated = true;
        }
      }
      //double gz = geomMan.GetGlobalZ( lnum );
      
      double tiltangle = sftclusterp->GetTiltAngle()*Deg2Rad;
      

      z.push_back( gz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( sftclusterp->GetLocalX() );
      ct.push_back( cos(tiltangle) ); st.push_back( sin(tiltangle) );

#if 0
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << sftclusterp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<sftclusterp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<sftclusterp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<sftclusterp->GetLocalHitPos()<< " " 
		<< std::endl;
#endif
    }
  }
  /*
  if(!(hitpattern[0] == true &&
      hitpattern[3] == true &&
      hitpattern[8] == true &&
      hitpattern[11] == true)) return status_ = false;
  */
  //nn = number of layers
  std::size_t nn = z.size();
  ///std::cout << "nn = " << nn << std::endl;

  double matrx[16], *mtp[4], fitp[4];
  mtp[0]=&matrx[0]; mtp[1]=&matrx[4]; mtp[2]=&matrx[8]; mtp[3]=&matrx[12];

  for( int i=0; i<4; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<4; ++j ){
      mtp[i][j]=0.0;
    }
  }
  
  //nn = number of layers
  for( std::size_t i=0; i<nn; ++i ){
    //weight, globalz , globalx, cos(tileangle) , sin(tiltangle)
    double ww=w[i], gz=z[i], ss=s[i], ctt=ct[i], stt=st[i];
    mtp[0][0] += ww*ctt*ctt;//cost*cost
    mtp[0][1] += ww*gz*ctt*ctt;//z*cost*cost
    mtp[0][2] += ww*ctt*stt;//cost*sint
    mtp[0][3] += ww*gz*ctt*stt;//z*cost*sint

    mtp[1][1] += ww*gz*gz*ctt*ctt;//z*z*cost*cost
    mtp[1][2] += ww*gz*ctt*stt;//z*cost*sint
    mtp[1][3] += ww*gz*gz*ctt*stt;//z*z*cost*sint

    mtp[2][2] += ww*stt*stt;//sint*sint
    mtp[2][3] += ww*gz*stt*stt;//z*sint*sint
    mtp[3][3] += ww*gz*gz*stt*stt;//z*z*sint*sint

    fitp[0] += ww*ss*ctt;//dl*cost
    fitp[1] += ww*gz*ss*ctt;//dl*z*cost
    fitp[2] += ww*ss*stt;//dl*cint
    fitp[3] += ww*gz*ss*stt;//dl*z*sint
  }
  mtp[1][0]=mtp[0][1]; mtp[2][0]=mtp[0][2]; mtp[3][0]=mtp[0][3];
  mtp[2][1]=mtp[1][2]; mtp[3][1]=mtp[1][3]; mtp[3][2]=mtp[2][3];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"   A1=  "<<fitp[0]<<"     A2=  "<<fitp[1]<<"     A3=  "<<fitp[2]
	   <<"     A4=  "<<fitp[3]<<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]<<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"    A24="<<mtp[1][3]<<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"    A34="<<mtp[2][3]<<std::endl;
  std::cout<<"    A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"    A44="<<mtp[3][3]<<std::endl;
#endif

  double Org[4][4]={{0}},Red[4][4]={{0}},Final[4][4]={{0}};
  double Org_vec[4]={0}, Solution_vec[4]={0};
  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Org[l][m]=mtp[l][m];
    }
  }

  for (int i=0; i<4; i++){
    Org_vec[i]=fitp[i];
  }

  if( MathTools::GaussJordan(mtp,4,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    //std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  x0_=fitp[0]; y0_=fitp[2]; u0_=fitp[1]; v0_=fitp[3];
  u1_=0.; v1_=0.;

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"     A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2]
	   <<"     A14="<<mtp[0][3]<<std::endl;
  std::cout<<"     A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2]
	   <<"     A24="<<mtp[1][3]<<std::endl;
  std::cout<<"     A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2]
	   <<"     A34="<<mtp[2][3]<<std::endl;
  std::cout<<"     A41="<<mtp[3][0]<<"    A42="<<mtp[3][1]<<"    A43="<<mtp[3][2]
	   <<"     A44="<<mtp[3][3]<<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"     A1="<<fitp[0]<<"    A2="<<fitp[1]<<"    A3="<<fitp[2]
	   <<"     A4="<<fitp[3]<<std::endl;
#endif


  for(int l=0; l<4;l++){
    for(int m=0; m<4; m++){
      Red[l][m]=mtp[l][m];
    }
  }

  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      Solution_vec[i] +=Red[i][j]*Org_vec[j];
    }
  }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]<<"    A2="<<Solution_vec[1]<<"    A3="<<Solution_vec[2]
	   <<"    A4="<<Solution_vec[3]<<std::endl;
#endif

  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      for(int k=0; k<4; k++){
	Final[i][j] += Red[i][k]*Org[k][j];
      }
    }
  }
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
    }
  }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"       A11="<<std::setw(10)<<Final[0][0]<<"       A12="<<std::setw(10)<<Final[0][1]
	   <<"       A13="<<std::setw(10)<<Final[0][2]<<"       A14="<<std::setw(10)<<Final[0][3]
	   <<std::endl;
  std::cout<<"       A21="<<std::setw(10)<<Final[1][0]<<"       A22="<<std::setw(10)<<Final[1][1]
	   <<"       A23="<<std::setw(10)<<Final[1][2]<<"       A24="<<std::setw(10)<<Final[1][3]
	   <<std::endl;
  std::cout<<"       A31="<<std::setw(10)<<Final[2][0]<<"       A32="<<std::setw(10)<<Final[2][1]
	   <<"       A33="<<std::setw(10)<<Final[2][2]<<"       A34="<<std::setw(10)<<Final[2][3]
	   <<std::endl;
  std::cout<<"       A41="<<std::setw(10)<<Final[3][0]<<"       A42="<<std::setw(10)<<Final[3][1]
	   <<"       A43="<<std::setw(10)<<Final[3][2]<<"       A44="<<std::setw(10)<<Final[3][3]
	   <<std::endl;

#endif

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], gz=z[i];
    double scal=GetX(gz)*ct[i]+GetY(gz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);

#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0_+u0_*gz)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*gz)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
  }
  chisqr /= nn-4.;
  /*  
  if(chisqr<2){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
  chisqr_=chisqr;
  for( std::size_t i=0; i<nn; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      int lnum = sftclusterp->GetLayer();
      double gz = GeomMan::GetInstance().GetGlobalZ( lnum );
      /*  
          if(chisqr<2){
          std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
          std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
          }
      */
      sftclusterp->SetProjectedPosition( GetX(gz), GetY(gz) );
    }
  }

  // std::cout << "***********************************************************" << std::endl;

  return status_=true;
}

bool TrLocalTrack::DoFit2( void )
{
  const std::string funcname = "[TrLocalTrack::DoFit2()]";
  const GeomMan & geomMan=GeomMan::GetInstance();
  
  std::size_t n = sftclusterArray_.size();
  
  if(n < TrLocalMinNHits2 ) return status_ = false;
  
  std::vector <double> z, w, s, ct, st;
  z.reserve(n); w.reserve(n); s.reserve(n);
  ct.reserve(n); st.reserve(n);
  
  for( std::size_t i=0; i<n; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      int lnum = sftclusterp->GetLayer();
      double ww = geomMan.GetResolution( lnum );
      double gz = geomMan.GetGlobalZ( lnum );
      double aa = sftclusterp->GetTiltAngle()*Deg2Rad;

      z.push_back( gz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( sftclusterp->GetLocalX() );//local x of the layer
      ct.push_back( cos(aa) );
      st.push_back( sin(aa) );
#if 0
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << sftclusterp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<sftclusterp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<sftclusterp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<sftclusterp->GetLocalHitPos()<< " " 
		<< std::endl;
#endif
    }
  }

  std::size_t nn = z.size();
  //std::cout << "nn = " << nn << std::endl;

  double matrx[36], *mtp[6], fitp[6];
  mtp[0]=&matrx[0]; 
  mtp[1]=&matrx[6]; 
  mtp[2]=&matrx[12]; 
  mtp[3]=&matrx[18];
  mtp[4]=&matrx[24];
  mtp[5]=&matrx[30];

  for( int i=0; i<6; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<6; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], gz=z[i], ss=s[i], ctt=ct[i], stt=st[i];

    mtp[0][0] += ww*ctt*ctt;//cost*cost
    mtp[0][1] += ww*gz*ctt*ctt;//z*cost*cost
    mtp[0][2] += ww*gz*gz*ctt*ctt;//z*z*cost*cost
    mtp[0][3] += ww*ctt*stt;//cost*sint
    mtp[0][4] += ww*gz*ctt*stt;//z*cost*sint
    mtp[0][5] += ww*gz*gz*ctt*stt;//z*z*cost*sint

    mtp[1][0] += ww*gz*ctt*ctt;//z*cost*cost
    mtp[1][1] += ww*gz*gz*ctt*ctt;//z*z*cost*cost
    mtp[1][2] += ww*gz*gz*gz*ctt*ctt;//z*z*z*cost*cost
    mtp[1][3] += ww*gz*ctt*stt;//z*cost*sint
    mtp[1][4] += ww*gz*gz*ctt*stt;//z*z*cost*sint
    mtp[1][5] += ww*gz*gz*gz*ctt*stt;//z*z*z*cost*sint

    mtp[2][0] += ww*gz*gz*ctt*ctt;//z*z*cost*cost
    mtp[2][1] += ww*gz*gz*gz*ctt*ctt;//z*z*z*cost*cost
    mtp[2][2] += ww*gz*gz*gz*gz*ctt*ctt;//z*z*z*z*cost*cost
    mtp[2][3] += ww*gz*gz*ctt*stt;//z*z*cost*sint
    mtp[2][4] += ww*gz*gz*gz*ctt*stt;//z*z*z*cost*sint
    mtp[2][5] += ww*gz*gz*gz*gz*ctt*stt;//z*z*z*z*cost*sint

    mtp[3][0] += ww*ctt*stt;//cost*sint
    mtp[3][1] += ww*gz*ctt*stt;//z*cost*sint
    mtp[3][2] += ww*gz*gz*ctt*stt;//z*z*cost*sint
    mtp[3][3] += ww*stt*stt;//sint*sint
    mtp[3][4] += ww*gz*stt*stt;//z*sint*sint
    mtp[3][5] += ww*gz*gz*stt*stt;//z*z*sint*sint

    mtp[4][0] += ww*gz*ctt*stt;//z*cost*sint
    mtp[4][1] += ww*gz*gz*ctt*stt;//z*z*cost*sint
    mtp[4][2] += ww*gz*gz*gz*ctt*stt;//z*z*z*cost*sint
    mtp[4][3] += ww*gz*stt*stt;//z*sint*sint
    mtp[4][4] += ww*gz*gz*stt*stt;//z*z*sint*sint
    mtp[4][5] += ww*gz*gz*gz*stt*stt;//z*z*z*sint*sint

    mtp[5][0] += ww*gz*gz*ctt*stt;//z*z*cost*sint
    mtp[5][1] += ww*gz*gz*gz*ctt*stt;//z*z*z*cost*sint
    mtp[5][2] += ww*gz*gz*gz*gz*ctt*stt;//z*z*z*z*cost*sint
    mtp[5][3] += ww*gz*gz*stt*stt;//z*z*sint*sint
    mtp[5][4] += ww*gz*gz*gz*stt*stt;//z*z*z*sint*sint
    mtp[5][5] += ww*gz*gz*gz*gz*stt*stt;//z*z*z*z*sint*sint

    fitp[0] += ww*ss*ctt;//dl*cost = x0_
    fitp[1] += ww*ss*gz*ctt;//dl*z*cost = u0_
    fitp[2] += ww*ss*gz*gz*ctt;//dl*z*z*cost = u1_

    fitp[3] += ww*ss*stt;//dl*sint = y0_
    fitp[4] += ww*ss*gz*stt;//dl*z*sint = v0_
    fitp[5] += ww*ss*gz*gz*stt;//dl*z*z*sint = v1_
  }

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"     A1=  "<<fitp[0]
	   <<"     A2=  "<<fitp[1]
	   <<"     A3=  "<<fitp[2]
	   <<"     A4=  "<<fitp[3]
	   <<"     A5=  "<<fitp[4]
	   <<"     A6=  "<<fitp[5]
	   <<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]
	   <<"    A12="<<mtp[0][1]
	   <<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]
	   <<"    A15="<<mtp[0][4]
      	   <<"    A16="<<mtp[0][5] <<std::endl;
  std::cout<<"    A21="<<mtp[1][0]
	   <<"    A22="<<mtp[1][1]
	   <<"    A23="<<mtp[1][2]
      	   <<"    A24="<<mtp[1][3]
	   <<"    A25="<<mtp[1][4]
      	   <<"    A26="<<mtp[1][5] <<std::endl;
  std::cout<<"    A31="<<mtp[2][0]
	   <<"    A32="<<mtp[2][1]
	   <<"    A33="<<mtp[2][2]
      	   <<"    A34="<<mtp[2][3]
	   <<"    A35="<<mtp[2][4]
      	   <<"    A36="<<mtp[2][5] <<std::endl;
  std::cout<<"    A41="<<mtp[3][0]
	   <<"    A42="<<mtp[3][1]
	   <<"    A43="<<mtp[3][2]
      	   <<"    A44="<<mtp[3][3]
	   <<"    A45="<<mtp[3][4]
      	   <<"    A46="<<mtp[3][5] <<std::endl;
  std::cout<<"    A51="<<mtp[4][0]
	   <<"    A52="<<mtp[4][1]
	   <<"    A53="<<mtp[4][2]
      	   <<"    A54="<<mtp[4][3]
	   <<"    A55="<<mtp[4][4]
      	   <<"    A56="<<mtp[4][5] <<std::endl;
  std::cout<<"    A61="<<mtp[5][0]
	   <<"    A62="<<mtp[5][1]
	   <<"    A63="<<mtp[5][2]
      	   <<"    A64="<<mtp[5][3]
	   <<"    A65="<<mtp[5][4]
      	   <<"    A66="<<mtp[5][5] <<std::endl;
#endif

  double Org[6][6]={{0}},Red[6][6]={{0}},Final[6][6]={{0}};
  double Org_vec[6]={0}, Solution_vec[6]={0};
  for(int l=0; l<6;l++){
    for(int m=0; m<6; m++){
      Org[l][m]=mtp[l][m];
    }
  }

  for (int i=0; i<6; i++){
    Org_vec[i]=fitp[i];
  }

  if( MathTools::GaussJordan(mtp,6,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  x0_=fitp[0]; u0_=fitp[1]; u1_=fitp[2];
  y0_=fitp[3]; v0_=fitp[4]; v1_=fitp[5];

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]
	   <<"    A12="<<mtp[0][1]
	   <<"    A13="<<mtp[0][2]
      	   <<"    A14="<<mtp[0][3]
	   <<"    A15="<<mtp[0][4]
      	   <<"    A16="<<mtp[0][5] <<std::endl;
  std::cout<<"    A21="<<mtp[1][0]
	   <<"    A22="<<mtp[1][1]
	   <<"    A23="<<mtp[1][2]
      	   <<"    A24="<<mtp[1][3]
	   <<"    A25="<<mtp[1][4]
      	   <<"    A26="<<mtp[1][5] <<std::endl;
  std::cout<<"    A31="<<mtp[2][0]
	   <<"    A32="<<mtp[2][1]
	   <<"    A33="<<mtp[2][2]
      	   <<"    A34="<<mtp[2][3]
	   <<"    A35="<<mtp[2][4]
      	   <<"    A36="<<mtp[2][5] <<std::endl;
  std::cout<<"    A41="<<mtp[3][0]
	   <<"    A42="<<mtp[3][1]
	   <<"    A43="<<mtp[3][2]
      	   <<"    A44="<<mtp[3][3]
	   <<"    A45="<<mtp[3][4]
      	   <<"    A46="<<mtp[3][5] <<std::endl;
  std::cout<<"    A51="<<mtp[4][0]
	   <<"    A52="<<mtp[4][1]
	   <<"    A53="<<mtp[4][2]
      	   <<"    A54="<<mtp[4][3]
	   <<"    A55="<<mtp[4][4]
      	   <<"    A56="<<mtp[4][5] <<std::endl;
  std::cout<<"    A61="<<mtp[5][0]
	   <<"    A62="<<mtp[5][1]
	   <<"    A63="<<mtp[5][2]
      	   <<"    A64="<<mtp[5][3]
	   <<"    A65="<<mtp[5][4]
      	   <<"    A66="<<mtp[5][5] <<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"    A1="<<fitp[0]
	   <<"    A2="<<fitp[1]
	   <<"    A3="<<fitp[2]
	   <<"    A4="<<fitp[3]
	   <<"    A5="<<fitp[4]
	   <<"    A6="<<fitp[5] <<std::endl;
#endif


  for(int l=0; l<6;l++){
    for(int m=0; m<6; m++){
      Red[l][m]=mtp[l][m];
    }
  }

  for(int i=0;i<6;i++){
    for(int j=0;j<6;j++){
      Solution_vec[i] +=Red[i][j]*Org_vec[j];
    }
  }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]
	   <<"    A2="<<Solution_vec[1]
	   <<"    A3="<<Solution_vec[2]
	   <<"    A4="<<Solution_vec[3]
	   <<"    A5="<<Solution_vec[4]
	   <<"    A6="<<Solution_vec[5]<<std::endl;
#endif

  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      for(int k=0; k<6; k++){
	Final[i][j] += Red[i][k]*Org[k][j];
      }
    }
  }
  for(int i=0; i<6; i++){
    for(int j=0; j<6; j++){
      if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
    }
  }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"       A11="<<std::setw(10)<<Final[0][0]
	   <<"       A12="<<std::setw(10)<<Final[0][1]
	   <<"       A13="<<std::setw(10)<<Final[0][2]
	   <<"       A14="<<std::setw(10)<<Final[0][3]
	   <<"       A15="<<std::setw(10)<<Final[0][4]
	   <<"       A16="<<std::setw(10)<<Final[0][5]
	   <<std::endl;
  std::cout<<"       A21="<<std::setw(10)<<Final[1][0]
	   <<"       A22="<<std::setw(10)<<Final[1][1]
	   <<"       A23="<<std::setw(10)<<Final[1][2]
	   <<"       A24="<<std::setw(10)<<Final[1][3]
	   <<"       A25="<<std::setw(10)<<Final[1][4]
	   <<"       A26="<<std::setw(10)<<Final[1][5]
	   <<std::endl;
  std::cout<<"       A31="<<std::setw(10)<<Final[2][0]
	   <<"       A32="<<std::setw(10)<<Final[2][1]
	   <<"       A33="<<std::setw(10)<<Final[2][2]
	   <<"       A34="<<std::setw(10)<<Final[2][3]
	   <<"       A35="<<std::setw(10)<<Final[2][4]
	   <<"       A36="<<std::setw(10)<<Final[2][5]
	   <<std::endl;
  std::cout<<"       A41="<<std::setw(10)<<Final[3][0]
	   <<"       A42="<<std::setw(10)<<Final[3][1]
	   <<"       A43="<<std::setw(10)<<Final[3][2]
	   <<"       A44="<<std::setw(10)<<Final[3][3]
	   <<"       A45="<<std::setw(10)<<Final[3][4]
	   <<"       A46="<<std::setw(10)<<Final[3][5]
	   <<std::endl;
  std::cout<<"       A51="<<std::setw(10)<<Final[4][0]
	   <<"       A52="<<std::setw(10)<<Final[4][1]
	   <<"       A53="<<std::setw(10)<<Final[4][2]
	   <<"       A54="<<std::setw(10)<<Final[4][3]
	   <<"       A55="<<std::setw(10)<<Final[4][4]
	   <<"       A56="<<std::setw(10)<<Final[4][5]
	   <<std::endl;
  std::cout<<"       A61="<<std::setw(10)<<Final[5][0]
	   <<"       A62="<<std::setw(10)<<Final[5][1]
	   <<"       A63="<<std::setw(10)<<Final[5][2]
	   <<"       A64="<<std::setw(10)<<Final[5][3]
	   <<"       A65="<<std::setw(10)<<Final[5][4]
	   <<"       A66="<<std::setw(10)<<Final[5][5]
	   <<std::endl;
#endif

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], gz=z[i];
    double scal=GetX(gz)*ct[i]+GetY(gz)*st[i];
    chisqr += ww*(s[i]-scal)*(s[i]-scal);

#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0_+u0_*gz)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*gz)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
  }
  chisqr /= nn-6.;
  /*  
  if(chisqr<2){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
  chisqr_=chisqr;
  for( std::size_t i=0; i<nn; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      int lnum = sftclusterp->GetLayer();
      double gz = GeomMan::GetInstance().GetGlobalZ( lnum );
      /*  
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
	  }
      */
      sftclusterp->SetProjectedPosition( GetX(gz), GetY(gz) );
      
    }
  }

  // std::cout << "***********************************************************" << std::endl;

  return status_=true;
}


bool TrLocalTrack::DoFitVXU( void )
{
  const std::string funcname = "[TrLocalTrack::DoFitVXU()]";
  const GeomMan & geomMan=GeomMan::GetInstance();

  std::size_t n = sftclusterArray_.size();

  if(n < TrLocalMinNHitsVXU ) return status_ = false;

  //double a=0,b=0;// <-Add !!
  double w[n+1],z[n+1],x[n+1];

  for( std::size_t i=0; i<n; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      int lnum = sftclusterp->GetLayer();
      w[i] = geomMan.GetResolution( lnum );
      z[i] = geomMan.GetGlobalZ( lnum );
      x[i] = sftclusterp->GetLocalX();

#if 0
      std::cout << "" << std::endl;
      std::cout << "**********" << std::endl;
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << sftclusterp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<sftclusterp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<sftclusterp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<sftclusterp->GetLocalHitPos()<< " " 
		<< std::endl;
      std::cout << "**********" << std::endl;
      std::cout << "" << std::endl;
#endif
    }
  }

  double A=0, B=0, C=0, D=0, E=0;// <-Add!! 
  for(std::size_t i=0; i<n; ++i){

    A += z[i]/(w[i]*w[i]);
    B += 1/(w[i]*w[i]);
    C += x[i]/(w[i]*w[i]);
    D += z[i]*z[i]/(w[i]*w[i]);
    E += x[i]*z[i]/(w[i]*w[i]);
  }

  a_ = (E*B-C*A)/(D*B-A*A);
  b_ = (D*C-E*A)/(D*B-A*A);
  c_ = 0.0;

  double chisqr = 0.0;

  for(std::size_t i=0; i<n; ++i){
    chisqr += (x[i]-a_*z[i]-b_)*(x[i]-a_*z[i]-b_)/(w[i]*w[i]);
  }

  chisqr /= n-2.;
  chisqr_=chisqr;

  //std::cout << "chisqr = " << chisqr << std::endl;
  //std::cout << "end "  << std::endl;
  //std::cout << " "  << std::endl;

  for( std::size_t i=0; i<n; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      //int lnum = sftclusterp->GetLayer();
      //double gz = GeomMan::GetInstance().GetGlobalZ( lnum );
      /*  
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
	  }
      */
     // sftclusterp->SetLocalCalPosVXU( a_*gz+b_ );
     std::cout << __FILE__ << "  " << __LINE__ << "this function is not implemeted yet ..." << std::endl;
    }
  }

  return status_=true;
}

bool TrLocalTrack::DoFitVXU2( void )
{
  const std::string funcname = "[TrLocalTrack::DoFitVXU2()]";
  const GeomMan & geomMan=GeomMan::GetInstance();
  
  std::size_t n = sftclusterArray_.size();
  
  if(n < TrLocalMinNHitsVXU2 ) return status_ = false;
  
  std::vector <double> z, w, s;
  z.reserve(n); w.reserve(n); s.reserve(n);
  
  for( std::size_t i=0; i<n; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      int lnum = sftclusterp->GetLayer();
      double ww = geomMan.GetResolution( lnum );
      double gz = geomMan.GetGlobalZ( lnum );

      z.push_back( gz ); w.push_back( 1./(ww*ww) ); 
      s.push_back( sftclusterp->GetLocalX() );

#if 0
      std::cout << std::setw(10) << "layer = " << lnum 
		<< std::setw(10) << "wire  = " << sftclusterp->GetWire() << " "
		<< std::setw(20) << "WirePosition = "<<sftclusterp->GetWirePosition() << " "
		<< std::setw(20) << "DriftLength = "<<sftclusterp->GetDriftLength() << " "
		<< std::setw(20) << "hit position = "<<sftclusterp->GetLocalHitPos()<< " " 
		<< std::endl;
#endif
    }
  }

  std::size_t nn = z.size();
  ///std::cout << "nn = " << nn << std::endl;

  double matrx[9], *mtp[3], fitp[3];
  mtp[0]=&matrx[0]; 
  mtp[1]=&matrx[3]; 
  mtp[2]=&matrx[6];

  for( int i=0; i<3; ++i ){
    fitp[i]=0.0;
    for( int j=0; j<3; ++j ){
      mtp[i][j]=0.0;
    }
  }

  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], gz=z[i], ss=s[i];
    mtp[0][0] += ww;//1
    mtp[0][1] += ww*gz;//z
    mtp[0][2] += ww*gz*gz;//z*z
    mtp[1][2] += ww*gz*gz*gz;//z*z*z
    mtp[2][2] += ww*gz*gz*gz*gz;//z*z*z*z

    fitp[0] += ww*ss;//s
    fitp[1] += ww*gz*ss;//z*s
    fitp[2] += ww*gz*gz*ss;//z*z*s
  }
  mtp[1][0]=mtp[0][1];
  mtp[1][1]=mtp[0][2];
  mtp[2][0]=mtp[0][2];
  mtp[2][1]=mtp[1][2];

  std::vector<int> indxc(nn), indxd(nn), ipiv(nn);

#if 0
  std::cout<<"             Vector:  Q_i               "<<std::endl;
  std::cout<<"     A1=  "<<fitp[0]
	   <<"     A2=  "<<fitp[1]
	   <<"     A3=  "<<fitp[2] <<std::endl;

  std::cout<<"             original matrix: M_ij        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2] <<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2] <<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2] <<std::endl;
#endif

  double Org[3][3]={{0}},Red[3][3]={{0}},Final[3][3]={{0}};
  double Org_vec[3]={0}, Solution_vec[3]={0};
  for(int l=0; l<3;l++){
    for(int m=0; m<3; m++){
      Org[l][m]=mtp[l][m];
    }
  }

  for (int i=0; i<3; i++){
    Org_vec[i]=fitp[i];
  }

  if( MathTools::GaussJordan(mtp,3,fitp,&indxc[0],
			     &indxd[0],&ipiv[0])==false ){
    std::cerr << funcname << ": Fitting fails" << std::endl;
    return status_=false;
  }
  b_=fitp[0]; a_=fitp[1]; c_=fitp[2];
  //f(z) = b_ + a_*x +c_*x*x

#if 0
  std::cout<<"             reduced matrix        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2] <<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2] <<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2] <<std::endl;
#endif

#if 0
  std::cout<<"             Solution from Gauss-Jordan             "<<std::endl;
  std::cout<<"    A1="<<fitp[0]
	   <<"    A2="<<fitp[1]
	   <<"    A3="<<fitp[2] <<std::endl;
#endif


  for(int l=0; l<3;l++){
    for(int m=0; m<3; m++){
      Red[l][m]=mtp[l][m];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Solution_vec[i] +=Red[i][j]*Org_vec[j];
    }
  }

#if 0
  std::cout<<"             Solution from explicite calculation       "<<std::endl;
  std::cout<<"    A1="<<Solution_vec[0]
	   <<"    A2="<<Solution_vec[1]
	   <<"    A3="<<Solution_vec[2] <<std::endl;
#endif

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	Final[i][j] += Red[i][k]*Org[k][j];
      }
    }
  }
  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      if(Final[i][j]<1.0e-10) Final[i][j]=0.0;
    }
  }

#if 0
  std::cout<< ""<<std::endl;
  std::cout<<"            final matrix        "<<std::endl;
  std::cout<<"    A11="<<mtp[0][0]<<"    A12="<<mtp[0][1]<<"    A13="<<mtp[0][2] <<std::endl;
  std::cout<<"    A21="<<mtp[1][0]<<"    A22="<<mtp[1][1]<<"    A23="<<mtp[1][2] <<std::endl;
  std::cout<<"    A31="<<mtp[2][0]<<"    A32="<<mtp[2][1]<<"    A33="<<mtp[2][2] <<std::endl;
#endif

  double chisqr=0.0;
  for( std::size_t i=0; i<nn; ++i ){
    double ww=w[i], gz=z[i];
    double scal= b_ + a_*gz +c_*gz*gz;
    chisqr += ww*(s[i]-scal)*(s[i]-scal);
    
#if 0
    if(1){
      //      std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
      //     std::cout<<"x coordinate = "<<(x0_+u0_*gz)<<std::endl;
      //     std::cout<<"y coordinate = "<<(y0_+v0_*gz)<<std::endl;
      std::cout<<std::setw(10)<<"layer = "<<i<<
	std::setw(10)<<"scal = "<<scal<<
	std::setw(10)<<"sdata = "<<s[i]<<std::endl;
      std::cout<<std::setw(10)<<"Res = "<<s[i]-scal<<std::endl;
      std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
      std::cout<<std::setw(10)<<"chisqr = "<<chisqr<<std::endl;
    }
#endif
  }
  chisqr /= nn-3.;
  /*  
  if(chisqr<2){
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
    std::cout << "chisqr = " << chisqr << " nn-4 = " << nn-4 << std::endl;
    std::cout<<" !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "<<std::endl;
  }
  */
  chisqr_=chisqr;
  for( std::size_t i=0; i<nn; ++i ){
    SFTCluster *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ){
      //int lnum = sftclusterp->GetLayer();
      //double gz = GeomMan::GetInstance().GetGlobalZ( lnum );
      /*  
	  if(chisqr<2){
	  std::cout<<std::setw(10)<<"lnum = "<< lnum <<std::endl;
	  std::cout<<std::setw(10)<<"X = "<< GetX(gz)<<" Y = "<< GetY(gz)<<std::endl;
	  }
      */
      //sftclusterp->SetCalPosition( GetX(gz), GetY(gz) );
      //sftclusterp->SetLocalCalPosVXU( b_ + a_*gz +c_*gz*gz );
      std::cout << __FILE__ << "  " << __LINE__ << "this function is not implemeted yet ..." << std::endl;
    }
  }

  // std::cout << "***********************************************************" << std::endl;

  return status_=true;
}

/*
bool TrLocalTrack::ReCalc( bool applyRecursively )
{
  static const std::string funcname = "[TrLocalTrack::ReCalc]";

  std::size_t n = sftclusterArray_.size();
  for( std::size_t i=0; i<n; ++i ){
    TrLTrackHit *sftclusterp = sftclusterArray_[i];
    if( sftclusterp ) sftclusterp->ReCalc( applyRecursively );
  }
  
  bool ret=DoFit();
  if( !ret ){
    std::cerr << funcname << ": Recalculation fails" << std::endl;
  }
  return ret;
}*/

