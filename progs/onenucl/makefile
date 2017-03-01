#
# Makefile for Glauber calculations:	Wim Cosyn
#
#-------------------------------------------------------------------------
#

#		Version number
VERSION =	

#		Calc-path
CPATH 	=	../..

#		C compilation...
#  CC      = 	icpc
CC      = 	g++

#		Flags
#  OPTFLAG	=	-O2 -fPIC
# OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel 

STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG)  #$(DBGFLAG) $(GPFLAG)  #$(STATIC)

#		Executable file
EFILE   = 	../../bin/JlabPAC16_calc ../../bin/test
EFILEST   = 	../../bin/staticonenucl


#		Mathematical functions
LIBS =	 	
WIMLIBS =   -lDIS -lKnockout -lGlauber -lMePhys -lIntegrator -lcuba
ROOT =		`root-config --cflags --glibs` -lMinuit2 `gsl-config --libs --cflags`
 LIBDIRS = -L $(CPATH)/lib
 INCLUDEDIRS =  -I $(CPATH)/Integrator -I $(CPATH)/Tools  -I $(CPATH)/Glauber -I $(CPATH)/DIS -I $(CPATH)/Knockout
STATICLIBS =  $(CPATH)/lib/libGlauber.a  $(CPATH)/lib/libDIS.a   $(CPATH)/lib/libMePhys.a $(CPATH)/lib/libKnockout.a  $(CPATH)/lib/libIntegrator.a


all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f *.o $(EFILE)

#-------------------------------------------------------------------------
.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp $(HDRS); 
	@echo $(@F).cpp 
	$(CC) -o $@ $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS)  $(LIBDIRS) $(ROOT)  $(INCLUDEDIRS) 
