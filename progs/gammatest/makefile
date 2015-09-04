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
  OPTFLAG	=	-O2 -fPIC
#OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel 

STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) # $(DBGFLAG) #$(GPFLAG)  #$(STATIC)

#		Executable file
EFILE   = 	../../bin/gammatest


#		Mathematical functions
LIBS =	 	-lm  
WIMLIBS = -lMePhys 
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/trunk/lib
 INCLUDEDIRS =  -I ~/Code/trunk/Integrator -I ~/Code/trunk/Tools  
STATICLIBS =   ~/Code/trunk/lib/libMePhys.a   ~/Code/trunk/lib/libIntegrator.a   


#		Header files
HDRS =		

# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o)  $(CPATH)/gamma.o

#-------------------------------------------------------------------------

#		Targets...

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) core*

#-------------------------------------------------------------------------

$(EFILE):	$(OBJS) $(HDRS)
		@echo "Linking $(EFILE)..."
		$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
$(EFILEST):	$(OBJS) $(HDRS) 
		@echo "Linking $(EFILEST)..."
		$(CC) $(CFLAGS) -o $(EFILEST) $(OBJS) $(STATICLIBS) $(LIBS) $(ROOT) $(LIBDIRS)


# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)

