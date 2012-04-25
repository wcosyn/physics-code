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
 OPTFLAG	=	-O3 -fPIC
#  OPTFLAG	=	-O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel 
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG)  #$(DBGFLAG)  $(GPFLAG) #$(STATIC)

#		Executable file
EFILE   = 	../../bin/rhoD
EFILEST   = 	../../bin/rhoDstatic

#		Mathematical functions
LIBS =	 	-lm  -lifcore
WIMLIBS = -lKnockout -lGlauber -lMePhys
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L/home/wim/Code/lib
 INCLUDEDIRS =  -I/home/wim/Code/Tools  -I/home/wim/Code/Glauber -I/home/wim/Code/Knockout
STATICLIBS =  /home/wim/Code/lib/libKnockout.a /home/wim/Code/lib/libGlauber.a  /home/wim/Code/lib/libMePhys.a 


#		Header files
HDRS =		

# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o)  $(CPATH)/rhod.o

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

