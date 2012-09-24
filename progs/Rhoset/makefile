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
  CC      = 	mpicxx
#  CC      = 	g++

#		Flags
#  OPTFLAG	=	-O3 -fPIC
 OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) # $(DBGFLAG)  $(GPFLAG) #$(STATIC)

#		Executable file
EFILE   = 	../../bin/rhofix
EFILEST   = 	../../bin/rhofixstatic

#		Mathematical functions
LIBS =	 	-lm  
WIMLIBS = -lKnockout -lGlauber -lMePhys -lIntegrator
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools  -I ~/Code/Glauber -I ~/Code/Knockout 
STATICLIBS =   ~/Code/lib/libKnockout.a  ~/Code/lib/libGlauber.a   ~/Code/lib/libMePhys.a ~/Code/lib/libIntegrator.a


#		Header files
HDRS =		

# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o)  $(CPATH)/rhoset.o

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

