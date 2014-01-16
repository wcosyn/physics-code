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
# CC      = 	icpc
 CC      = 	g++

#		Flags
 OPTFLAG	=	-g #-O3 -fPIC
#OPTFLAG	=	 -O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG)  #$(DBGFLAG)  #$(GPFLAG) #$(STATIC)


#		Mathematical functions
LIBS =	 	-lm  #-lifcore
WIMLIBS = -lMePhys -lGlauber -lDIS -lIntegrator -lKnockout -lcuba
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L$(HOME)/Code/lib
 INCLUDEDIRS =  -I$(HOME)/Code/Tools -I$(HOME)/Code/Glauber -I$(HOME)/Code/DIS -I$(HOME)/Code/Integrator -I$(HOME)/Code/Knockout
STATICLIBS = $(HOME)/Code/lib/libglauber.a  $(HOME)/Code/lib/libmephys.a $(HOME)/Code/lib/libdis.a $(HOME)/Code/lib/libcuba.a
LDFLAGS = -Wl,-rpath,$(HOME)/Code/lib/
#		Header files
HDRS =		

# root headers

#		Object files
#OBJS = 	        $(HDRS:.hpp=.o) 

#-------------------------------------------------------------------------

#		Targets...

targets=cubatest cubatest_wrapped

all:		$(targets)

$(targets): %: %.cpp
	$(CC) $(CFLAGS) $< $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(LIBDIRS) $(INCLUDEDIRS) -o $@
	
#-------------------------------------------------------------------------
.PHONY:clean
clean:;		@rm $(targets) *.o

#-------------------------------------------------------------------------



#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
