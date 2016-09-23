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
OPTFLAG	=	-O3 -fPIC -std=c++11
# OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel 
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) #$(DBGFLAG) #$(GPFLAG)    #$(STATIC)

#		Executable file
EFILE   = 	../../bin/deepstest ../../bin/deepsfit ../../bin/bonustest ../../bin/bonusfit  ../../bin/bonusextrapolate ../../bin/bonusextrapolate2 ../../bin/bonus_F2nfit ../../bin/bonus_F2nfit_allangles ../../bin/bonusfitreadout ../../bin/bonusfit_rig ../../bin/bonus_F2nfit_sys ../../bin/b1_calc ../../bin/hermes_comp ../../bin/b1_calc_comp
EFILEST   = 	../../bin/staticdeepstest ../../bin/staticdeepsfit ../../bin/staticbonustest ../../bin/staticbonusfit  ../../bin/staticbonusscatterfit


#		Mathematical functions
LIBS =	 	
WIMLIBS =   -lDIS -lKnockout -lGlauber -lMePhys -lIntegrator -lcuba
ROOT =		`root-config --cflags --glibs` -lMinuit2 `gsl-config --libs --cflags`
 LIBDIRS = -L $(CPATH)/lib
 INCLUDEDIRS =  -I $(CPATH)/Integrator -I $(CPATH)/Tools  -I $(CPATH)/Glauber -I $(CPATH)/DIS -I $(CPATH)/Knockout
STATICLIBS =  $(CPATH)/lib/libGlauber.a  $(CPATH)/lib/libDIS.a   $(CPATH)/lib/libMePhys.a $(CPATH)/lib/libKnockout.a  $(CPATH)/lib/libIntegrator.a

#		Header files
HDRS =		bonusdata.h bonusfits.h bonusfits2.h

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f *.o $(EFILE)

#-------------------------------------------------------------------------
.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp $(HDRS); 
	@echo $(@F).cpp 
	$(CC) -o $@ $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS)  $(LIBDIRS) $(ROOT)  $(INCLUDEDIRS) 

