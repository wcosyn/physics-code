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
OPTFLAG	=	-O3 -mtune=core2 -march=corei7-avx -fp-model source -ipo -fPIC
STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) # $(DBGFLAG)  $(GPFLAG) #$(STATIC)

#		Executable file
EFILE   = 	$(CPATH)/../lib/libdis.so.1.0
EFILEST   = 	$(CPATH)/../lib/libdis.a


#		Mathematical functions
LIBS =	 	-lm -lifcore
WIMLIBS = -lmephys -lglauber
ROOT =		`/usr/local/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L/home/wim/Code/lib
 INCLUDEDIRS =  -I/home/wim/Code/Tools  -I/home/wim/Code/Glauber
STATICLIBS = /home/wim/Code/lib/libglauber.a  /home/wim/Code/lib/libmephys.a 

#		Header files
HDRS =		$(CPATH)/InclusiveCross.hpp\
		$(CPATH)/DeuteronMomDistr.hpp\
		$(CPATH)/DeuteronCross.hpp


# root headers

#		Object files
OBJS = 	        $(HDRS:.hpp=.o)

#-------------------------------------------------------------------------

#		Targets...

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f $(OBJS) core*

#-------------------------------------------------------------------------

$(EFILE):	$(OBJS) $(HDRS) $(DICTHDRS)
		@echo "Linking library $(EFILE)..."
		$(CC) -shared -Wl,-soname,libdis.so.1 -o $(EFILE) $(OBJS) $(LIBS) $(ROOT) $(LIBDIRS)
		ln -sf ../lib/libdis.so.1.0 ../lib/libdis.so
		ln -sf ../lib/libdis.so.1.0 ../lib/libdis.so.1
    
$(EFILEST):	$(OBJS) $(HDRS) $(DICTHDRS) $(HDRS2) $(INCHDRS)
		@echo "Linking static library $(EFILEST)..."
		xiar -cvr $(EFILEST) $(OBJS)


# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
%.o: 		%.cpp 
		@echo "Compiling..."
		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)
