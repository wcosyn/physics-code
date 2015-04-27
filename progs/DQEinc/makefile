#
# Makefile for Glauber calculations:	Wim Cosyn
#
#-------------------------------------------------------------------------
#

#		Version number
VERSION =	

#		Calc-path
CPATH 	=	.

#   CC      = 	icpc
 CC      = 	g++

#		Flags
 OPTFLAG	=	-O3 -fPIC -std=c++0x
#   OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -std=c++0x
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) $(DBGFLAG)  #$(GPFLAG) #$(STATIC)

#		Executable file
EFILE1   = 	../../bin/DQEinc


#		Mathematical functions
LIBS =	 	-lm 
WIMLIBS = -lDIS -lKnockout -lGlauber -lMePhys -lIntegrator -lcuba
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/trunk/lib
 INCLUDEDIRS =  -I ~/Code/trunk/Integrator -I ~/Code/trunk/Tools  -I ~/Code/trunk/Glauber -I ~/Code/trunk/Knockout -I ~/Code/trunk/DIS
STATICLIBS =   ~/Code/trunk/lib/libDIS.a ~/Code/trunk/lib/libKnockout.a  ~/Code/trunk/lib/libGlauber.a   ~/Code/trunk/lib/libMePhys.a ~/Code/trunk/lib/libIntegrator.a

#		Header files
HDRS =		

# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o) $(CPATH)/DQEinc.o

#-------------------------------------------------------------------------

#		Targets...

all:		$(EFILE1)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) core*

#-------------------------------------------------------------------------

$(EFILE1):	$(OBJS) $(HDRS)
		@echo "Linking $(EFILE)..."
		$(CC) $(CFLAGS) -o $(EFILE1) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)

# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
