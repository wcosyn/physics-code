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
# CC      = 	g++

#		Flags
# OPTFLAG	=	-O3 -fPIC
OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel 
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) #$(GPFLAG)  #$(DBGFLAG)  #$(STATIC)

#		Executable file
EFILE   = 	../../bin/deepsfit
EFILEST   = 	../../bin/staticdeepsfit


#		Mathematical functions
LIBS =	 	-lm  -lifcore
WIMLIBS = -lMePhys -lGlauber -lDIS -lKnockout -lIntegrator
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L/home/wim/Code/lib
 INCLUDEDIRS =  -I/home/wim/Code/Integrator -I/home/wim/Code/Tools  -I/home/wim/Code/Glauber -I/home/wim/Code/DIS -I/home/wim/Code/Knockout
STATICLIBS = /home/wim/Code/lib/libglauber.a  /home/wim/Code/lib/libmephys.a /home/wim/Code/lib/libdis.a /home/wim/Code/lib/libknockout.a /home/wim/Code/lib/libintegrator.a

#		Header files
HDRS =		
# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o) $(CPATH)/deepsfit.o

#-------------------------------------------------------------------------

#		Targets...

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) core*

#-------------------------------------------------------------------------

$(EFILE):	$(OBJS) $(HDRS)
		@echo "Linking $(EFILE)..."
		$(CC) $(CFLAGS) -o $(EFILE) $(MAIN1) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
$(EFILEST):	$(OBJS) $(HDRS) 
		@echo "Linking $(EFILEST)..."
		$(CC) $(CFLAGS) -o $(EFILEST) $(MAIN1) $(OBJS) $(STATICLIBS) $(LIBS) $(ROOT) $(LIBDIRS)

# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
