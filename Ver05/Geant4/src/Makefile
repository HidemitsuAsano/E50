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
EVTGENDIR	= /home/sirotori/Generators/EvtGen/R01-01-00
EVTGENINCDIR	= $(EVTGENDIR)
EVTGENLIBDIR	= $(EVTGENDIR)/lib
EVTGENINCLUDE	= -I$(EVTGENINCDIR)
EVTGENLIB	= -L$(EVTGENLIBDIR) -lEvtGen -lEvtGenExternal
EVTGEN_EXTERNAL = 0
EXTRAFLAGS      = -DEVTGEN_EXTERNAL
#
#
HEPMCDIR	= /usr/local
HEPMCINCDIR	= $(HEPMCDIR)/include/HepMC
HEPMCLIBDIR	= $(HEPMCDIR)/lib
HEPMCINCLUDE	= -I$(HEPMCINCDIR)
HEPMCLIB	= -L$(HEPMCLIBDIR) -lHepMC
#
#
PYTHIADIR	= /home/sirotori/Generators/PYTHIA/pythia8170
PYTHIAINCDIR	= $(PYTHIADIR)/include
PYTHIALIBDIR	= $(PYTHIADIR)/lib
PYTHIAINCLUDE	= -I$(PYTHIAINCDIR)
PYTHIALIB	= -L$(PYTHIALIBDIR) -lpythia8 -llhapdfdummy -lhepmcinterface
#
#
PHOTOSDIR	= /home/sirotori/Generators/PHOTOS
PHOTOSINCDIR	= $(PHOTOSDIR)/include
PHOTOSLIBDIR	= $(PHOTOSDIR)/lib
PHOTOSINCLUDE	= -I$(PHOTOSINCDIR)
PHOTOSLIB	= -L$(PHOTOSLIBDIR) -lPhotosCxxInterface -lPhotosFortran
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
			$(EVTGENINCLUDE) $(HEPMCINCLUDE)		\
			$(PYTHIAINCLUDE) $(PHOTOSINCLUDE)		\
			$(BOOSTINCLUDE)
#
LIBS		= $(G4LIBS) $(CLHEPLIB) $(MESALIBS) $(XLIBS)		\
			$(CERNLIBS) $(RTLIBS)	$(EVTGENLIB) $(HEPMCLIB)	\
			$(PYTHIALIB) $(PHOTOSLIB) $(BOOSTLIB)
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
		SpecField.o				\
		SpecFieldMap.o				\
		SimpleFieldElement.o			\
		MaterialList.o				\
		DetectorConstruction.o			\
		s_T0Wall.o                              \
                s_BFTU.o                                \
                s_BFTD.o                                \
                s_SFT.o                                 \
                s_AFT.o                                 \
                s_ITWC1.o                               \
                s_ITWC2L.o                              \
                s_ITWC2R.o                              \
                s_STWC1.o                               \
                s_STWC2.o                               \
                s_TofWall.o                             \
                s_ITofWall.o                            \
                s_ITofSWall.o                           \
                s_PAD.o                                 \
                s_RICH.o                                \
                s_PID1.o                                \
                s_PID2.o                                \
                s_MuonFilter.o                          \
		r_T0Wall.o                              \
		TargetHit.o				\
		TargetSD.o				\
		T0Hit.o	         			\
		T0SD.o	               			\
		BTrackerHit.o				\
		BTrackerSD.o	               		\
		SITrackerHit.o				\
		SITrackerSD.o	               		\
		AOTrackerHit.o				\
		AOTrackerSD.o	               		\
		IITrackerHit.o				\
		IITrackerSD.o	               		\
		IOTrackerHit.o	  			\
		IOTrackerSD.o	               		\
		SOTrackerHit.o				\
		SOTrackerSD.o	               		\
		TofHit.o				\
		TofSD.o	               			\
		ITofHit.o				\
		ITofSD.o	               		\
		PADHit.o				\
		PADSD.o 	               		\
		RICHHit.o				\
		RICHSD.o 	               		\
		PID1Hit.o				\
		PID1SD.o 	               		\
		PID2Hit.o				\
		PID2SD.o 	               		\
		MuonFilterHit.o				\
		MuonFilterSD.o 	               		\
		VDHit.o	         			\
		VDSD.o	               			\
		PhysicsList.o				\
		RunAction.o				\
		PrimaryGeneratorAction.o		\
		GeneratorHelper.o			\
		AngDisGenerator.o			\
		PGReaction1.o				\
		PGReaction2.o				\
		PGReaction3.o				\
		PGReaction4.o				\
		PGReaction5.o				\
		PGReaction6.o				\
		Legendre.o				\
		EventAction.o				\
		ReadJAMData.o				\
		JAMData.o				\
		SteppingAction.o			\
		Analysis.o				\
		Particle.o				\
		VisManager.o				\
		Transportation.o			\
		PidUtil.o				\
		LiqTarget.o

OBJS_SPEC =	MainSpec.o				\
		DetectorConstructionSpec.o		\
		AnalysisSpec.o			        \
		EventActionSpec.o			


all: $(BINDIR)/simSpec 

clean:
	$(RM) $(OBJS) $(OBJS_SPEC) *.d *~
#
SRCS	= $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

$(BINDIR)/simSpec: $(OBJS_SPEC) $(OBJS)
	$(C++) -o $@ $(CFLAGS) $(OBJS_SPEC) $(OBJS) $(LIBS)

