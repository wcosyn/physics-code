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
 OPTFLAG	=	-O3 -fPIC -std=c++11
# OPTFLAG	=	-E -O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel 

STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG)  $(DBGFLAG) #$(GPFLAG)  #$(STATIC)
HOME    = /home/wim

#		Mathematical functions
LIBS =	 	-lm
WIMLIBS = -lMePhys -lGlauber -lKnockout -lIntegrator 
ROOT =		`root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L$(HOME)/Code/trunk/lib
 INCLUDEDIRS =  -I$(HOME)/Code/trunk/Integrator -I$(HOME)/Code/trunk/Tools  -I$(HOME)/Code/trunk/Glauber -I $(HOME)/Code/trunk/Knockout
STATICLIBS =  $(HOME)/Code/trunk/lib/libKnockout.a  $(HOME)/Code/trunk/lib/libGlauber.a  $(HOME)/Code/trunk/lib/libMePhys.a   $(HOME)/Code/trunk/lib/libIntegrator.a   


#		Header files
HDRS =		

# root headers

#		Object files
# OBJS = 	        $(EFILE:../../bin/=)  
# OBJS2 = 	        $(HDRS:.hpp=.o)  $(CPATH)/neutestint.o

#-------------------------------------------------------------------------
#		Executable file
EFILE   = 	../../bin/neutest ../../bin/neutestint ../../bin/neutestfixEin
EFILEST   = 	../../bin/staticneutest 

#		Targets...

all:		$(EFILE)
static:		$(EFILEST)
#-------------------------------------------------------------------------
clean:;		@rm -f *.o $(EFILE)

#-------------------------------------------------------------------------
.SECONDEXPANSION:
$(EFILE):  $$(@F).cpp; 
	@echo $(@F).cpp
	$(CC) $(CFLAGS) $(@F).cpp $(LIBDIRS) $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(INCLUDEDIRS) -o $@

# $(EFILE):; @echo $(@F).cpp

# # $(EFILE):	$(OBJS) $(HDRS)
# # 		@echo "Linking $(EFILE)..."
# # 		$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
# $(EFILE2):	$(OBJS2) $(HDRS)
# 		@echo "Linking $(EFILE2)..."
# 		$(CC) $(CFLAGS) -o $(EFILE) $(OBJS) $(LIBS) $(WIMLIBS) $(ROOT) $(LIBDIRS)
# $(EFILEST):	$(OBJS) $(HDRS) 
# 		@echo "Linking $(EFILEST)..."
# 		$(CC) $(CFLAGS) -o $(EFILEST) $(OBJS) $(STATICLIBS) $(LIBS) $(ROOT) $(LIBDIRS)


# $(OBJS):	$(HDRS) $(HDRS2) $(DICTHDRS) $(INCHDRS)

#------


#		Rule to make objectfile
# %.o: 		%.cpp
# 		@echo "Compiling..."
# 		$(CC) $(CFLAGS) -c $*.cpp -o $*.o $(ROOT) $(INCLUDEDIRS)

