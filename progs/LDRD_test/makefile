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
EFILE   = 	../../bin/LDRD_test 
EFILEST   = 	../../bin/staticLDRD_test 


#		Mathematical functions
LIBS =	 	
WIMLIBS =   -lDIS -lKnockout -lMePhys -lIntegrator -lcuba
ROOT =		`root-config --cflags --glibs` -lMinuit2  `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools -I ~/Code/DIS -I ~/Code/Knockout
STATICLIBS =   ~/Code/lib/libmephys.a  ~/Code/lib/libdis.a  ~/Code/lib/libknockout.a  ~/Code/lib/libintegrator.a

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

