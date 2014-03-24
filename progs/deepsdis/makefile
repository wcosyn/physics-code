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
#  CC      = 	icpc
CC      = 	g++

#		Flags
OPTFLAG	=	-O3 -fPIC
# OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel 
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) #$(GPFLAG)  #$(DBGFLAG)  #$(STATIC)

#		Executable file
EFILE   = 	../../bin/deepstest ../../bin/deepsfit ../../bin/bonustest
EFILEST   = 	../../bin/staticdeepstest ../../bin/staticdeepsfit ../../bin/staticbonustest


#		Mathematical functions
LIBS =	 	-lm 
WIMLIBS = -lMePhys -lGlauber -lDIS -lKnockout -lIntegrator -lcuba
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools  -I ~/Code/Glauber -I ~/Code/DIS -I ~/Code/Knockout
STATICLIBS =  ~/Code/lib/libglauber.a   ~/Code/lib/libmephys.a  ~/Code/lib/libdis.a  ~/Code/lib/libknockout.a  ~/Code/lib/libintegrator.a

#		Header files
HDRS =		

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f *.o $(EFILE)

#-------------------------------------------------------------------------
.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp; 
	@echo $(@F).cpp
	$(CC) $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(LIBDIRS) $(INCLUDEDIRS) -o $@

