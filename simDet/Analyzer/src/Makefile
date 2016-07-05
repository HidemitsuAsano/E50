##	Makefile
#
#
### CLHEP
CLHEPDIR        = /usr/local
#
CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR     = $(CLHEPDIR)/lib
#
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIBS       = -L$(CLHEPLIBDIR) -lCLHEP
#
#
### RooT
RTINCLUDE	= `root-config --cflags`
RTLIBS		= `root-config --libs`
#
#
CFLAGS      = -g -O2 -Wall
#-Wall 
#DFLAGS      = -Df2cFortran -Dextname -DDEBUG 

#debug_flags += -DMemoryLeak

FLAGS = $(CFLAGS) $(DFLAGS) -I. $(CLHEPINCLUDE) $(RTINCLUDE)

LIBOBJ = libSpecAna.a
SPECLIB = -L. -lSpecAna
LIBS = $(CLHEPLIBS) $(RTLIBS)
LIBSONL = $(CLHEPLIBS) $(RTLIBS) 

CC = gcc
C++ = g++

BINDIR = ../bin

SRCS    = $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

OBJ =   ConfMan.o RawData.o \
	GetNumberFromKernelEntropyPool.o \
        PrimInfo.o HodoRawHit.o RawHit.o \
        TrGeomMan.o TrGeomRecord.o \
	Kinematics.o MathTools.o \
	HodoAnalyzer.o HodoHit.o HodoCluster.o\
	TrHit.o TrAnalyzer.o TrTrackSearch.o \
	TrHitCluster.o TrLocalTrack.o TrLTrackHit.o 

lib: $(OBJ)
	ar r $(LIBOBJ) $(OBJ)

all: analysMonitor \
     analysBeamTracking \

.cc.o:  
	$(C++) $(FLAGS) -MMD -c $<

analysMonitor: lib UserMonitor.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserMonitor.o 	\
	Main.o $(SPECLIB) $(LIBS)

analysBeamTracking: lib UserBeamTracking.o Main.o $(LIBOBJ)
	$(C++) $(FLAGS) -o $(BINDIR)/$@ UserBeamTracking.o 	\
	Main.o $(SPECLIB) $(LIBS)

clean:
	rm -f *.o *.d *~ core* libSpecAna.a ../bin/*
