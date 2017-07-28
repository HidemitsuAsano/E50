#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <string>
#include <vector>
#include <signal.h>

#include "ConfMan.hh"
#include "VEvent.hh"

#include <TROOT.h>
#include <TFile.h>

enum eArg
  {
    kArgProcessName,
    kArgConfFile,
    kArgInFile,
    kArgOutRootFile,
    kArgc
  };

const int MaxChar = 400;

void closeFile(int sig)
{
  if(gFile)
    {
      gFile->Write();
      gFile->Close();
    }
}


TROOT theROOT("sample-ana", "test");

int main(int argc, char* argv[])
{
  if (argc!=kArgc)
    {
      std::cout <<" #D Usage: "
		<<" [analyzer config file]"
		<<" [Input data (.dat file)] "
		<<" [output root file]"
		<< std::endl;
      return 0;
    }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  std::ifstream InputData;
  InputData.open(argv[kArgInFile],std::ios::in|std::ios::binary);
  const std::string& rootFile = arg[kArgOutRootFile];
  
  
  //open pedestal table
  std::string sdata = argv[kArgInFile];
  int ipath = sdata.find_last_of("/")+1;
  int iext  = sdata.find_last_of(".");
  std::string tablename = sdata.substr(ipath,iext-ipath);
  std::string table_dir = "pede_table/";
  tablename +="_ret.txt";
  std::cout << "table name: " << tablename.c_str() << std::endl;
  table_dir += tablename;
  //std::ifstream ftable;
  //ftable.open(table_dir.c_str());
  //if(!ftable){
   // std::cout << "NOT FOUND : " <<  table_dir.c_str() << std::endl;
  //}


  ConfMan* gconfManager =new ConfMan(confFile,table_dir);

  signal(SIGINT,closeFile);

  TFile f(rootFile.c_str(), "recreate");
  std::cout << "#D recreate root file : "
	    << rootFile << std::endl;

  if (!gconfManager->Initialize())
    return 0;
  std::cout << "ConfMan::AnaMode " << gconfManager->AnaMode() << std::endl;
  if( !InputData ) return 0;
  int evNum=0;
  
  

  while( InputData ){
   // VEvent* event = gconfManager->EventAllocator();
    VEvent* event = gconfManager->EventAllocator();
    if ( event->ProcessingNormal( InputData ) ){
      delete event;
    }else{ 
      std::cout << "Analysis finished ! " << std::endl;
      break;
    }
    ++evNum;
    if( evNum%500 == 0 ){
      std::cout << "EventNum=" << evNum << std::endl;
    }
  }

  gFile->Write();
  gFile->Close();
  
  return 0;
}

