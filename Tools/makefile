#
# Makefile for Glauber calculations:	Wim Cosyn
#
#-------------------------------------------------------------------------
#

#		Version number
VERSION =	

#		Calc-path
CPATH 	=	.

#		C compilation...
 CC      = 	icpc
 FF	=	ifort
# CC      = 	g++
# FF	=	gfortran

#		Flags
# OPTFLAG	=	-O3 -fPIC
 OPTFLAG	=	-O3 -mtune=core2 -march=corei7-avx -fp-model source -ipo -unroll -pthread -fPIC
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG)   #$(DBGFLAG)  #$(GPFLAG) #$(STATIC)

#		Executable file
EFILE   = 	$(CPATH)/../lib/libmephys.so.1.0
EFILEST   = 	$(CPATH)/../lib/libmephys.a


#		Mathematical functions
LIBS =	 	-lm 
ROOT =		`/usr/local/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
# LIBDIRS = -L/home/wim/OOPGlauber/src
# INCLUDEDIRS = -I/home/wim/OOPGlauber/src


#		Header files
HDRS =		$(CPATH)/NucleonEMOperator.hpp\
		$(CPATH)/NuclStructure.hpp\
		$(CPATH)/DeuteronStructure.hpp\
		$(CPATH)/Utilfunctions.hpp
		
HDRS2=		$(CPATH)/TCstPWF.h\
		$(CPATH)/TDeuteron.h\
		$(CPATH)/TGrossPWF.h\
		$(CPATH)/TInterpolatingWavefunction.h\
		$(CPATH)/TWavefunctionImplementation.h\
		$(CPATH)/TYukawaPWF.h\
		$(CPATH)/Matrix.h\
		$(CPATH)/GammaStructure.h\
		$(CPATH)/FourVector.h\
		$(CPATH)/TSpinor.h\
		$(CPATH)/TLorentzQuaternion.h\
		$(CPATH)/TKinematics.h\
		$(CPATH)/TKinematics2to2.h\
		$(CPATH)/TElectronKinematics.h\
		$(CPATH)/pair.h\
		$(CPATH)/recmosh.h\
		$(CPATH)/newcoef.h\
		$(CPATH)/TMPI.h\
		$(CPATH)/TKinematics2to3.h\
		$(CPATH)/TKinematics2to3WithLabAngles.h

INCHDRS = 	$(CPATH)/constants.hpp

# root headers
DICTHDRS = 	$(CPATH)/TCstPWFDict.h\
		$(CPATH)/TDeuteronDict.h\
		$(CPATH)/TGrossPWFDict.h\
		$(CPATH)/TInterpolatingWavefunctionDict.h\
		$(CPATH)/TWavefunctionImplementationDict.h\
		$(CPATH)/TYukawaPWFDict.h\
		$(CPATH)/TSpinorDict.h\
		$(CPATH)/TLorentzQuaternionDict.h\
		$(CPATH)/TElectronKinematicsDict.h\
		$(CPATH)/TKinematicsDict.h\
		$(CPATH)/TKinematics2to2Dict.h\
		$(CPATH)/TMPIDict.h\
		$(CPATH)/TKinematics2to3Dict.h\
		$(CPATH)/TKinematics2to3WithLabAnglesDict.h


FORTFILES = 	$(CPATH)/alekhin.f\
		$(CPATH)/F1F209_test.f

#		Object files
OBJS = 	        $(HDRS:.hpp=.o) $(HDRS2:.h=.o) $(DICTHDRS:.h=.o) $(FORTFILES:.f=.o)


#-------------------------------------------------------------------------

#		Targets...

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) $(DICTHDRS) $(DICTHDRS:.h=.cpp) core*

#-------------------------------------------------------------------------

$(EFILE):	$(OBJS) $(HDRS) $(DICTHDRS) $(HDRS2) $(INCHDRS)
		@echo "Linking library $(EFILE)..."
		$(CC) -shared -Wl,-soname,libmephys.so.1 -o $(EFILE) $(OBJS) $(LIBS) $(ROOT)
		ln -sf ../lib/libmephys.so.1.0 ../lib/libmephys.so
		ln -sf ../lib/libmephys.so.1.0 ../lib/libmephys.so.1
$(EFILEST):	$(OBJS) $(HDRS) $(DICTHDRS) $(HDRS2) $(INCHDRS)
		@echo "Linking static library $(EFILEST)..."
		xiar -cvr $(EFILEST) $(OBJS)
# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT)
%Dict.cpp:	%.h %LinkDef.h
		@echo "Making dictionary..."
		echo "@rootcint -f $*Dict.cpp -c -p `root-config --cflags` $(CFLAGS) $*.h $*LinkDef.h"
		@rootcint -f $*Dict.cpp -c -p $*.h $*LinkDef.h

$(FORTFILES:.f=.o): 	$(FORTFILES)
		@echo "Compiling..."
		$(FF) -fPIC -c $*.f
