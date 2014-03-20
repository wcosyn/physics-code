#
# Makefile for Glauber calculations:	Wim Cosyn
#
#-------------------------------------------------------------------------
#

#		Version number
VERSION =	

#		Calc-path
CPATH 	=	.

#CC      = 	icpc
 CC      = 	g++

#		Flags
 OPTFLAG	=	-O3 -fPIC
#  OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) $(DBGFLAG)  #$(GPFLAG) #$(STATIC)

#		Executable file
# EFILE1   = 	../../bin/inclusive
# EFILE2   = 	../../bin/disinc
# EFILE3   = 	../../bin/disuniform
# EFILE4   = 	../../bin/f2n
# EFILEST   = 	../../bin/staticinclusive


#		Mathematical functions
LIBS =	 	-lm 
WIMLIBS = -lDIS -lKnockout -lGlauber -lMePhys -lIntegrator -lcuba
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools  -I ~/Code/Glauber -I ~/Code/Knockout -I ~/Code/DIS
STATICLIBS =  ~/Code/lib/libDIS.a ~/Code/lib/libKnockout.a  ~/Code/lib/libGlauber.a  ~/Code/lib/libMePhys.a   ~/Code/lib/libIntegrator.a   

#		Header files
HDRS =		

# root headers

#		Object files
#		Object files
# OBJS = 	        $(EFILE:../../bin/=)  

#-------------------------------------------------------------------------
#		Executable file
EFILE   = 	../../bin/Azzres ../../bin/Azzuniform ../../bin/inclusive

#		Targets...

all:		$(EFILE)
clean:;		@rm -f *.o $(EFILE)

.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp; 
	@echo $(@F).cpp
	$(CC) $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(LIBDIRS) $(INCLUDEDIRS) -o $@
#-------------------------------------------------------------------------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
