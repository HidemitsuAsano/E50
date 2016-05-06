#
#	Makefile
#
# Compilers
C++	= g++
CC	= gcc
#
# geant4.9.3
GEANT4DIR	= /cern/geant4
#
G4LIBDIR	= $(GEANT4DIR)/lib/Linux-g++
G4INCDIR	= $(GEANT4DIR)/include
#
G4INCLUDE	= -I$(G4INCDIR)
G4LIBS		= -L$(G4LIBDIR)			\
			-lG4global		\
			-lG4digits_hits		\
			-lG4geometry		\
			-lG4interfaces		\
			-lG4parmodels		\
			-lG4particles		\
			-lG4processes		\
			-lG4event		\
			-lG4graphics_reps	\
			-lG4intercoms		\
			-lG4materials		\
			-lG4modeling		\
			-lG4readout		\
			-lG4run			\
			-lG4track		\
			-lG4tracking		\
			-lG4vis_management	\
			-lG4FR			\
			-lG4VRML		\
			-lG4OpenGL		\
			-lG4RayTracer		\
			-lG4visXXX		\
			-lG4gl2ps		
#                       -lG4emstandard          
#
CLHEPDIR	= /usr/local
CLHEPINCDIR	= $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR	= $(CLHEPDIR)/lib
CLHEPINCLUDE	= -I$(CLHEPINCDIR)
CLHEPLIB	= -L$(CLHEPLIBDIR) -lCLHEP
#
#
X11DIR		= /usr/X11R6
X11LIBDIR	= $(X11DIR)/lib
X11INCDIR	= $(X11DIR)/include
#
XINCLUDE	= -I$(X11INCDIR)
XLIBS		= -L$(X11LIBDIR) -lX11 -lXmu
#
MESADIR		= /usr/X11R6/lib
MESALIBS	= -L$(MESADIR) -lGL -lGLU
#
#
BOOSTDIR	= /usr
BOOSTLIBDIR	= $(BOOSTDIR)/lib64
BOOSTINCDIR	= $(BOOSTDIR)/include
BOOSTINCLUDE	= -I$(BOOSTINCDIR)
BOOSTLIB	= -L$(BOOSTLIBDIR) -lboost_iostreams
#
#
RTINCLUDE	= `root-config --cflags`
RTLIBS		= `root-config --libs`
#
DEFINEFLAGS	=	-DG4UI_USE_TCSH		\
			-DG4UI_USE_TERMINAL	\
			-DG4VIS_USE		\
			-DG4VIS_USE_OPENGL	\
			-DG4VIS_USE_OPENGLX	\
			-DG4VIS_USE_VRMLFILE	\
			-DG4VIS_USE_DAWNFILE	\
			-DG4VIS_USE_RAYTRACER	\
			-DG4USE_STD_NAMESPACE	\
			-DEVTGEN_EXTERNAL
#
DEBUGFLAGS	= -g -Wall
#
CPPFLAGS	= $(DEFINEFLAGS) $(DEBUGFLAGS) -I. $(G4INCLUDE)		\
			$(CLHEPINCLUDE) $(XINCLUDE) $(RTINCLUDE)	\
			$(BOOSTINCLUDE)
#
LIBS		= $(G4LIBS) $(CLHEPLIB) $(MESALIBS) $(XLIBS)		\
			$(CERNLIBS) $(RTLIBS) $(BOOSTLIB)
#
CFLAGS		= -O2
#
FLAGS		= $(CFLAGS) $(CPPFLAGS)
#
BINDIR		= ../bin
#
TOUCH		= touch
#
RM		= rm -f
#
# Rules
#
.cc.o:
	$(C++) $(FLAGS) -MMD -c $<
.cpp.o:
	$(C++) $(FLAGS) -MMD -c $<
#
.c.o:
	$(CC) $(FLAGS) -MMD -c $<
#
#

OBJS	= 	GetNumberFromKernelEntropyPool.o	\
		ConfMan.o				\
		DCGeomMan.o				\
		DCGeomRecord.o				\
		MaterialList.o				\
		DetectorConstruction.o			\
                s_SFT.o                                 \
                s_T0Wall.o                              \
                r_SFT_A.o                               \
                r_SFT_B.o                               \
                r_SFT_C.o                               \
		SFTHit.o				\
		SFTSD.o 	               		\
		T0Hit.o		         		\
		T0SD.o	                   		\
		PhysicsList.o				\
		RunAction.o				\
		PrimaryGeneratorAction.o		\
		GeneratorHelper.o			\
		EventAction.o				\
		Transportation.o 			\
		SteppingAction.o			\
		Analysis.o				\
		VisManager.o				\
		PidUtil.o				\

OBJS_SPEC =	MainSpec.o				\
		DetectorConstructionSpec.o		\
		AnalysisSpec.o			        \
		EventActionSpec.o			


all: $(BINDIR)/simDet

clean:
	$(RM) $(OBJS) $(OBJS_SPEC) *.d *~ $(BINDIR)/simDet
#
SRCS	= $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

$(BINDIR)/simDet: $(OBJS_SPEC) $(OBJS)
	$(C++) -o $@ $(CFLAGS) $(OBJS_SPEC) $(OBJS) $(LIBS)

