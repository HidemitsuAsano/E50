/*
  Main.cc

  2012/5  K.Shirotori
*/

#include <iostream>
#include <fstream>
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
		<<" [data input stream] "
		<<" [output root file]"
		<< std::endl;
      return 0;
    }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  std::ifstream InputData(argv[kArgInFile]);
  const std::string& rootFile = arg[kArgOutRootFile];

  ConfMan* gconfManager =new ConfMan( confFile );

  signal(SIGINT,closeFile);

  TFile f(rootFile.c_str(), "recreate");
  std::cout << "#D recreate root file : "
	    << rootFile << std::endl;

  if (!gconfManager->Initialize())
    return 0;

  if( !InputData ) return 0;
  int evNum=0;
  while( !InputData.eof() ){
    VEvent* event = gconfManager->EventAllocator();
    if ( event->ProcessingNormal( InputData ) )
      delete event;
    else break;
    ++evNum;
    if( evNum%1000 == 0 ){
      std::cout << "EventNum=" << evNum << std::endl;
    }
  }

  gFile->Write();
  gFile->Close();
  
  return 0;
}

