TARGET1  = anaMissMass01
#TARGET2 = anaMissMass02
#TARGET3 = anaMissMass03
#TARGET4 = anaMissMass10
#TARGET5 = anaMissMass20
#TARGET6 = anaMissMass21
#TARGET7 = anaMissMass30

all: $(TARGET1) #$(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) $(TARGET7)

ROOTFLAG += -DUSE_ROOT

ROOT_OBJ1  = anaMissMass01.o Particle.o
#ROOT_OBJ2 = anaMissMass02.o Particle.o
#ROOT_OBJ3 = anaMissMass03.o Particle.o
#ROOT_OBJ4 = anaMissMass10.o Particle.o
#ROOT_OBJ5 = anaMissMass20.o Particle.o
#ROOT_OBJ6 = anaMissMass21.o Particle.o
#ROOT_OBJ7 = anaMissMass30.o Particle.o

ROOTSYS  = /usr/local/cern/root
CFLAGS = -Dlynx
CC     = /usr/bin/g++
LD     = /usr/bin/gcc
EXTRALIBS =  -L/usr/X11R6/lib -lX11 -lnsl -lcrypt -ldl -lc -lm 
MACHINE = LINUX

### CLHEP
CLHEPDIR        = /usr/local
CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR     = $(CLHEPDIR)/lib
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIBS       = -L$(CLHEPLIBDIR) -lCLHEP

anaMissMass01.o: anaMissMass01.cc Makefile
#anaMissMass02.o: anaMissMass02.cc Makefile
#anaMissMass03.o: anaMissMass03.cc Makefile
#anaMissMass10.o: anaMissMass10.cc Makefile
#anaMissMass20.o: anaMissMass20.cc Makefile
#anaMissMass21.o: anaMissMass21.cc Makefile
#anaMissMass30.o: anaMissMass30.cc Makefile

ifeq ($(ROOTFLAG),-DUSE_ROOT)
# use root-config command
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
CFLAGS  += $(ROOTCFLAGS) $(ROOTFLAG) -I. $(CLHEPINCLUDE)
LIBS    += $(ROOTLIBS) $(CLHEPLIBS)
OBJ1 += $(ROOT_OBJ1)
OBJ2 += $(ROOT_OBJ2)
OBJ3 += $(ROOT_OBJ3)
OBJ4 += $(ROOT_OBJ4)
OBJ5 += $(ROOT_OBJ5)
OBJ6 += $(ROOT_OBJ6)
OBJ7 += $(ROOT_OBJ7)
endif


CFLAGS += -D$(MACHINE) -DOFFLINE $(OPT)
CFLAGS += -O2
CFLAGS += -g
CFLAGS += -I./include

LIBS += $(EXTRALIBS)

$(TARGET1): $(OBJ1)
	$(LD) -o $(TARGET1) $(OBJ1) $(LIBS)

#$(TARGET2): $(OBJ2)
#	$(LD) -o $(TARGET2) $(OBJ2) $(LIBS)

#$(TARGET3): $(OBJ3)
#	$(LD) -o $(TARGET3) $(OBJ3) $(LIBS)

#$(TARGET4): $(OBJ4)
#	$(LD) -o $(TARGET4) $(OBJ4) $(LIBS)

#$(TARGET5): $(OBJ5)
#	$(LD) -o $(TARGET5) $(OBJ5) $(LIBS)

#$(TARGET6): $(OBJ6)
#	$(LD) -o $(TARGET6) $(OBJ6) $(LIBS)

#$(TARGET7): $(OBJ7)
#	$(LD) -o $(TARGET7) $(OBJ7) $(LIBS)

.SUFFIXES: .cc

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@
.cc.o:
	$(CC) $(CFLAGS) -c $< -o $@
clean:
	rm -f $(TARGET1) core *.o 
#	rm -f $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) $(TARGET6) core *.o 
