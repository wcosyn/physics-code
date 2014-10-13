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
PROJECTPATH 	= 	../..

#		Flags
OPTFLAG	=	-O3 -fPIC
# OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel 
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) #$(GPFLAG)  #$(DBGFLAG)  #$(STATIC)

#		Executable file
EFILE   = 	../../bin/LDRD_test 
EFILEST   = 	../../bin/staticLDRD_test 


#		Mathematical functions
LIBS =	 	
WIMLIBS =   -lDIS -lKnockout -lMePhys -lIntegrator -lcuba
ROOT =		`root-config --cflags --glibs` -lMinuit2  `gsl-config --cflags --libs`  
 LIBDIRS = -L $(PROJECTPATH)/lib
 INCLUDEDIRS =  -I $(PROJECTPATH)/Integrator -I $(PROJECTPATH)/Tools -I $(PROJECTPATH)/DIS -I $(PROJECTPATH)/Knockout
STATICLIBS =   $(PROJECTPATH)/lib/libmephys.a  $(PROJECTPATH)/lib/libdis.a  $(PROJECTPATH)/lib/libknockout.a  $(PROJECTPATH)/lib/libintegrator.a

#		Header files
# HDRS =		bonusdata.h bonusfits.h

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f *.o $(EFILE)

#-------------------------------------------------------------------------
.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp $(HDRS); 
	@echo $(@F).cpp 
	$(CC) -o $@ $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(LIBDIRS) $(INCLUDEDIRS) 

