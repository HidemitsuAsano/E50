#include "GeomMan.hh"
#include "GeomRecord.hh"
#include "DetectorID.hh"

#include <string>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <cstdlib>

// chpattern_: member variable to change channel map
// 0: default ch. map
//    0-63
//
// 1: correct flat cable pin assign
//   0-15 -> 15-0
//  16-31 -> 31-16
//  32-47 -> 47-32
//  48-63 -> 63-48
//
// 2: rotate MPPC board
//   0-31 -> 31-0
//  32-63 -> 63-32
//
// 3: rotate MPPC board (another roration)
//   0-15 -> 16-31
//  16-31 ->  0-15
//  32-47 -> 48-63
//  48-63 -> 32-47


const int MaxChar = 200;

GeomMan *GeomMan::geomMan_=0;

GeomMan::GeomMan():
chpattern_(1)
{
   if(chpattern_ > 0){
     std::cout << "GeomMan:: " <<  "ch. map pattern : "<< chpattern_ << std::endl;
   }
}

GeomMan::~GeomMan()
{}


GeomMan & GeomMan::GetInstance( void )
{
  if( !geomMan_ ){
    geomMan_ = new GeomMan();
  }
  return *geomMan_;
}

double GeomMan::GetGlobalZ( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetGlobalZ(int)]"; 
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->length_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetResolution( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetResolution(int)]"; 
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->resol_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetTiltAngle( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetTiltAngle(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->tiltAngle_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetRotAngle1( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetRotAngle1(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->rotAngle1_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetRotAngle2( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetRotAngle2(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->rotAngle2_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

double GeomMan::GetOffset(int lnum) const 
{
  static const std::string funcname = "[GeomMan::GetOffset(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->offset_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

/*
const ThreeVector & GeomMan::GetGlobalPosition( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetGlobalPosition(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->pos_;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector GeomMan::NormalVector( int lnum ) const
{
  static const std::string funcname = "[GeomMan::NormalVector(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->NormalVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}

ThreeVector GeomMan::UnitVector( int lnum ) const
{
  static const std::string funcname = "[GeomMan::UnitVector(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo->UnitVector();
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}*/

const GeomRecord * GeomMan::GetRecord( int lnum ) const
{
  static const std::string funcname = "[GeomMan::GetRecord(int)]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ) return pGeo;
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 

//calculate ch ( = fiber) local x position from geant geometry. 
//GeomRecord::ChPos(int ch) is modified from original code.
//here, "ch" means the fiber ID
double GeomMan::calcChPosition( int lnum, double ch ) const
{
  static const std::string funcname = "[GeomMan::calcChPosition()]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ){
    return pGeo->ChPos(ch);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}


int GeomMan::calcChNumber( int lnum, double pos ) const
{
  static const std::string funcname = "[GeomMan::calcChNumber()]";
  int id = GetSFTID(lnum);
  GeomRecord *pGeo = geomRecord_[id];
  if( pGeo ){
    return pGeo->ChNumber(pos);
  }
  else{
    std::cerr << funcname << ": No record. Layer#=" 
	      << lnum << std::endl;
    throw std::out_of_range(funcname+": No record" );
  }
}
 

void GeomMan::clearElements( void )
{
  //  for_each( geomRecord_.begin(), geomRecord_.end(), DeleteObject() );
  std::map <int, GeomRecord *>::iterator itr;
  for( itr=geomRecord_.begin(); itr!=geomRecord_.end(); ++itr ){
    delete itr->second;
  }
  geomRecord_.clear();
}


bool GeomMan::Initialize( void )
{
  static const std::string funcname = "[GeomMan::Initialize]";
  char str[MaxChar];
  char cname[MaxChar];
  int id;
  double xs, ys, zs, ta, ra1, ra2, l, res, w0, dd, ofs;

  FILE *fp;

  if( ( fp = fopen( filename_.c_str(), "r" ) ) == 0 ){
    throw std::invalid_argument(funcname+": file open fail");
  }

  clearElements();

  while( fgets( str, MaxChar, fp ) != 0 ){
    if( str[0]!='#' ){
      if( sscanf( str, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		  &id, cname, &xs, &ys, &zs, &ta, &ra1, &ra2, &l, &res,
		  &w0, &dd, &ofs )
	  == 13 ){
	GeomRecord *pRec = 
	  new GeomRecord( id, cname, xs, ys, zs, 
			    ta, ra1, ra2, l, res,
			    w0, dd, ofs );

	GeomRecord *pOld = geomRecord_[id];
	geomRecord_[id] = pRec;

	if( pOld ){
	  std::cerr << funcname << ": duplicated id number. "
		    << " following record is deleted." << std::endl;
	  std::cerr << "Id=" << pOld->id_ << " " <<  " ) ... " << std::endl;
	  delete pOld;
	}
      }
      else {
	std::string strtemp=str;
	std::cerr << funcname << ": Invalid format " << strtemp << std::endl;
      }
    }
  }

  fclose(fp);

  return true;
}

int GeomMan::GetSFTID(int lnum = -1) const
{
  //range check
  if(lnum <PlMinSFT || lnum > PlMaxSFT){
    std::cerr << __FILE__ << "  " << __LINE__ << " invalid SFT layer number! " << lnum  << std::endl;
    return -1;
  }
  int sftID = lnum + 1;

  return sftID;
}


std::vector <int> GeomMan::GetDetectorIDList( void ) const
{
  std::vector<int> vlist;
  vlist.reserve(geomRecord_.size());
  std::map <int, GeomRecord *>::const_iterator 
    itr=geomRecord_.begin(), end=geomRecord_.end();

  for(; itr!=end; ++itr ){
    vlist.push_back( itr->first );
  }

  return vlist;
}

//proto type BFT ch. scheme 
//0 origin, up to 127 ch.
//MPPC board 1: ch  0 - 31  L0 subL0
//MPPC board 2: ch 32 - 63  L1 subL1
//MPPC board 3: ch 64 - 95  L2 subL0
//MPPC board 4: ch 96 - 127 L3 subL1
//layer 0 (beam downstream) - 3 (beam upstream)
//EASIROC 0-3
//
//fiber no. is defined in each layer
//0 origin, up to fiber 11.
//
//X:0,U:1,V:2
//layer configuration
//<-downstream    upstream->
//X,U,V,X,U,V,U,V,X,U,V,X


//return X:0,U:1,V:2
int GeomMan::getXUV(int ch = -1) const
{
  
  if(64<= ch) ch -= 64;
  
  ch = ChPattern(ch);

  if(0 <= ch && ch <= 7){
    return 1;   //U (8 fibers)
  }else if(8  <= ch && ch <= 17){
    return 1;  //U (10 fibers)
  }else if(18 <= ch && ch <= 23){
    return 0;  //X (6 fibers)
  }else if(24 <= ch && ch <= 25){
    return 1;  //U (2 fibers)
  }else if(26 <= ch && ch <= 31){
    return 0; //X (6 fibers)
  }else if(32 <= ch && ch <= 37){
    return 0; //X (6 fibers)
  }else if(38 <= ch && ch <= 39){
    return 2; //V (2 fibers)
  }else if(40 <= ch && ch <= 45){
    return 0; //X (6 fibers)
  }else if(46 <= ch && ch <= 63){
    return 2; //V (18 fibers)
  }else{ 
    std::cout << "invalid ch: " << ch << std::endl;
    return -1;
  }
}

//return sublayer  0 or 1
int GeomMan::getsublayer(int ch = -1) const
{
  if(64<= ch) ch -= 64;
  
  ch = ChPattern(ch);

  int type = getXUV(ch);
  int sublayer = -1;
  if(type == 1){//U
    if(0 <= ch && ch <= 7) sublayer = 0; //L0 U0-U8 
    else if(8  <= ch && ch <= 15) sublayer = 1;//L1 U0-U7
    else if(16 <= ch && ch <= 17) sublayer = 0;//L0 U8-U9
    else if(24 <= ch && ch <= 25) sublayer = 1;//L1 U8-U9
  }else if(type == 0){//X
    if(18 <= ch && ch <= 23) sublayer = 0; //L0 X0-X5
    else if(26 <= ch && ch <= 31) sublayer = 1;//L1 X0-X5
    else if(32 <= ch && ch <= 37) sublayer = 0;//L0 X6-X11
    else if(40 <= ch && ch <= 45) sublayer = 1;//L1 X6-X11                               
  }else if(type == 2){//V
    if(38 <= ch && ch <= 39) sublayer = 0 ; //L0 V0-V1
    else if(46 <= ch && ch <= 47) sublayer = 1 ;//L1 V0-V1
    else if(48 <= ch && ch <= 55) sublayer = 0 ;//L0 V2-V9 
    else if(56 <= ch && ch <= 63) sublayer = 1 ; //L1 V2-V9
  }else{
    std::cout << "invalid ch" << ch << std::endl;
  }
  
  return sublayer;
}

//return layer number 0-11
//<--downstream    upstream-->
//0  1  2  3  4  5  6  7  8  9 10 11
//X  U  V  X  U  V  U  V  X  U  V  X
int GeomMan::getlayer(int ch = -1) const
{
  bool is64over = false;
  if(64<= ch){
    ch -= 64;
    is64over = true;
  }
  
  int type = getXUV(ch);
  
  ch = ChPattern(ch);

  int layer = -1;
  if(type==1){//U
    if(0 <= ch && ch <= 7) layer = 1; //L0 U0-U8 
    else if(8 <= ch && ch <= 15) layer = 4;//L1 U0-U7
    else if(16<= ch && ch <= 17) layer = 1;//L0 U8-U9
    else if(24 <= ch && ch <= 25) layer = 4;//L1 U8-U9
    
    if(is64over) layer+=5;
  
  }else if(type == 0){//X
    if(18<= ch && ch <= 23) layer = 0; //L0 X0-X5
    else if(26<= ch && ch <= 31) layer = 3;//L1 X0-X5
    else if(32<= ch && ch <= 37) layer = 0;//L0 X6-X11
    else if(40<= ch && ch <= 45) layer = 3;//L1 X6-X11      
    
    if(is64over) layer+=8;

  }else if(type == 2){//V
    if(38<= ch && ch <= 39) layer = 2 ; //L0 V0-V1
    else if(46<= ch && ch <= 47) layer = 5 ;//L1 V0-V1
    else if(48<= ch && ch <= 55) layer = 2 ;//L0 V2-V9 
    else if(56<= ch && ch <= 63) layer = 5 ; //L1 V2-V9
    
    if(is64over) layer+=5;
  
  }else{
    std::cout << "invalid ch" << ch << std::endl;
  }

  return layer;
}

int GeomMan::getfiber(int ch = -1) const
{
  if(64<= ch) ch = ch - 64;
  int type = getXUV(ch);
  
  ch = ChPattern(ch);

  int fiber = -1;
  //U type
  if(type==1){//U
    if(0 <= ch && ch <= 7) fiber = ch; // L0 U0-U8 
    else if(8 <= ch && ch <= 15) fiber = ch - 8;//L1 U0-U7
    else if(16<= ch && ch <= 17) fiber = ch - 8;//L0 U8-U9
    else if(24 <= ch && ch <= 25) fiber = ch - 16;//L1 U8-U9
  }else if(type == 0){//X
    if(18<= ch && ch <= 23) fiber = ch - 18; //L0 X0-X5
    else if(26<= ch && ch <= 31) fiber = ch - 26;//L1 X0-X5
    else if(32<= ch && ch <= 37) fiber = ch - 26;//L0 X6-X11
    else if(40<= ch && ch <= 45) fiber = ch - 34;//L1 X6-X11                               
  }else if(type == 2){//V
    if(38<= ch && ch <= 39) fiber = ch - 38 ; //L0 V0-V1
    else if(46<= ch && ch <= 47) fiber = ch - 46;//L1 V0-V1
    else if(48<= ch && ch <= 55) fiber = ch - 46;//L0 V2-V9 
    else if(56<= ch && ch <= 63) fiber = ch - 54; //L1 V2-V9
  }else{
    std::cout << "invalid ch" << ch << std::endl;
  }
   
  return fiber;
}


//return layer0 -3
int GeomMan::getbiglayer(int ch = -1) const
{
  bool is64over = false;
  if(64<= ch){
    ch -= 64;
    is64over = true;
  }
  
  int type = getXUV(ch);
  
  ch =  ChPattern(ch);

  int layer = -1;
  if(type==1){//U
    if(0 <= ch && ch <= 7) layer = 0; //L0 U0-U8 
    else if(8 <= ch && ch <= 15) layer = 1;//L1 U0-U7
    else if(16<= ch && ch <= 17) layer = 0;//L0 U8-U9
    else if(24 <= ch && ch <= 25) layer = 1;//L1 U8-U9
  }else if(type == 0){//X
    if(18<= ch && ch <= 23) layer = 0; //L0 X0-X5
    else if(26<= ch && ch <= 31) layer = 1;//L1 X0-X5
    else if(32<= ch && ch <= 37) layer = 0;//L0 X6-X11
    else if(40<= ch && ch <= 45) layer = 1;//L1 X6-X11                               
  }else if(type == 2){//V
    if(38<= ch && ch <= 39) layer = 0 ; //L0 V0-V1
    else if(46<= ch && ch <= 47) layer = 1 ;//L1 V0-V1
    else if(48<= ch && ch <= 55) layer = 0 ;//L0 V2-V9 
    else if(56<= ch && ch <= 63) layer = 1 ; //L1 V2-V9
  }else{
    std::cout << "invalid ch" << ch << std::endl;
  }
  if(is64over) layer += 2; 
  return layer;
}



//get easiroc chip #0-3 
int GeomMan::geteasiroc(int ch = -1) const
{

  ch = ChPattern(ch);

  int easiroc = -1;
  if(0 <= ch && ch <= 31) easiroc = 0;
  else if(32<= ch && ch <= 63) easiroc = 1;
  else if(64<= ch && ch <= 95) easiroc = 2;
  else if(96<= ch && ch <= 127) easiroc = 3;
  else std::cout << "ivalid ch" << ch << std::endl;
  return easiroc;

}


int GeomMan::ChPattern(int ch = -1) const 
{

  if(chpattern_ == 0){
    return ch;
  }else if(chpattern_ == 1){
    //  0-15 -> 15-0
    // 16-31 -> 31-16
    // 32-47 -> 47-32
    // 48-63 -> 63-48
    if(0<= ch && ch <=15)  ch = 15 - ch; 
    else if(16<= ch && ch<=31) ch = 47 - ch;
    else if(32<= ch && ch<=47) ch = 79 - ch;
    else if(48<= ch && ch<=63) ch = 111- ch;
    else{
      std::cout << "Invalid channel number !: " << ch << std::endl;
    }
  }else if(chpattern_ == 2){
    //  0-31 -> 31-0
    // 32-63 -> 63-32
    if(0<=ch && ch<= 31) ch = 31 - ch;
    else if(32<= ch && ch <=63) ch = 95 - ch;
    else{
      std::cout << "Invalid channel number !: " << ch << std::endl;
    }
  }else if(chpattern_ == 3){
    //  0-15 -> 16-31
    // 16-31 ->  0-15
    // 32-47 -> 48-63
    // 48-63 -> 32-47
    if(0<=ch && ch<=15) ch  = ch + 16;
    else if(16<=ch && ch<=31) ch = ch - 16;
    else if(32<=ch && ch<=47) ch = ch + 16;
    else if(48<=ch && ch<=63) ch = ch - 16;
    else{
      std::cout << "Invalid channel number !: " << ch << std::endl;
    }
  }else{
    std::cout << "invalid ch. pattern !: " << chpattern_ << std::endl;
  }

  return ch;
}



