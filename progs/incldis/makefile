#
# Makefile for Glauber calculations:	Wim Cosyn
#
#-------------------------------------------------------------------------
#

#		Version number
VERSION =	

#		Calc-path
CPATH 	=	.

 CC      = 	icpc
#  CC      = 	g++

#		Flags
#  OPTFLAG	=	-O3 -fPIC
 OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) $(DBGFLAG)  #$(GPFLAG) #$(STATIC)

#		Executable file
EFILE1   = 	../../bin/inclusive
EFILE2   = 	../../bin/disinc
EFILE3   = 	../../bin/disuniform
EFILEST   = 	../../bin/staticinclusive


#		Mathematical functions
LIBS =	 	-lm 
WIMLIBS = -lDIS -lKnockout -lGlauber -lMePhys -lIntegrator
ROOT =		`$(ROOTSYS)/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools  -I ~/Code/Glauber -I ~/Code/Knockout -I ~/Code/DIS
STATICLIBS =   ~/Code/lib/libDIS.a ~/Code/lib/libKnockout.a  ~/Code/lib/libGlauber.a   ~/Code/lib/libMePhys.a ~/Code/lib/libIntegrator.a

#		Header files
HDRS =		

# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o) $(CPATH)/inclusive.o
OBJS2 = 	        $(HDRS:.hpp=.o) $(CPATH)/disinctest.o
OBJS3 = 	        $(HDRS:.hpp=.o) $(CPATH)/disuniform.o

#-------------------------------------------------------------------------

#		Targets...

# all:		$(EFILE1)
inc:	$(EFILE1)
disinc:	$(EFILE2)
uniform: $(EFILE3)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) core*

#-------------------------------------------------------------------------

$(EFILE1):	$(OBJS) $(HDRS)
		@echo "Linking $(EFILE)..."
		$(CC) $(CFLAGS) -o $(EFILE1) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
$(EFILE2):	$(OBJS2) $(HDRS)
		@echo "Linking $(EFILE2)..."
		$(CC) $(CFLAGS) -o $(EFILE2) $(OBJS2) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
$(EFILE3):	$(OBJS3) $(HDRS)
		@echo "Linking $(EFILE3)..."
		$(CC) $(CFLAGS) -o $(EFILE3) $(OBJS3) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
$(EFILEST):	$(OBJS) $(HDRS) 
		@echo "Linking $(EFILEST)..."
		$(CC) $(CFLAGS) -o $(EFILEST) $(OBJS) $(STATICLIBS) $(LIBS) $(ROOT) $(LIBDIRS)

# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
