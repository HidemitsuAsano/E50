##	Makefile
#
#
### CLHEP
#CLHEPDIR        = /usr/local
#
#CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
#CLHEPLIBDIR     = $(CLHEPDIR)/lib
#
#CLHEPINCLUDE    = -I$(CLHEPINCDIR)
#CLHEPLIBS       = -L$(CLHEPLIBDIR) -lCLHEP
#
#
### RooT
RTINCLUDE	= `root-config --cflags`
RTLIBS		= `root-config --libs`
#
#
CFLAGS      = -g -O2 -Wall -Wextra
#-Wall 
#DFLAGS      = -Df2cFortran -Dextname -DDEBUG 

#debug_flags += -DMemoryLeak

#FLAGS = $(CFLAGS) $(DFLAGS) -I. $(CLHEPINCLUDE) $(RTINCLUDE)
FLAGS = $(CFLAGS) $(DFLAGS) -I. $(RTINCLUDE)

LIBOBJ = libSpecAna.a
SPECLIB = -L. -lSpecAna
#LIBS = $(CLHEPLIBS) $(RTLIBS)
LIBS = $(RTLIBS)
#LIBSONL = $(CLHEPLIBS) $(RTLIBS) 
LIBSONL = $(RTLIBS) 

CC = gcc
C++ = g++
#C++ = clang

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o RawData.o \
	GetNumberFromKernelEntropyPool.o \
        SFTRawHit.o \
        TrGeomMan.o TrGeomRecord.o \
	MathTools.o \
	TrHit.o TrAnalyzer.o TrTrackSearch.o \
	SFTCluster.o TrLocalTrack.o  

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: anaBFT 

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

anaBFT: lib AnaBFT.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ AnaBFT.o 	\
	Main.o $(SPECLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libSpecAna.a ../bin/*
