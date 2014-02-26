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
 OPTFLAG	=	-O3 -fPIC
# OPTFLAG	=	-E -O3 -fPIC -shared-intel -fp-model source -unroll -pthread -ipo -xsse3 -parallel 

STATIC  =	--static
DBGFLAG	=	-g
GPFLAG  = 	-pg
WFLAG   = 	-Wall
CFLAGS	=	$(OPTFLAG) $(WFLAG) # $(DBGFLAG) #$(GPFLAG)  #$(STATIC)


#		Mathematical functions
LIBS =	 	-lm
WIMLIBS = -lMePhys -lGlauber -lKnockout -lIntegrator -lcuba
ROOT =		`$(ROOTSYS)/bin/root-config --cflags --glibs` -lMathMore `gsl-config --cflags --libs`  
 LIBDIRS = -L ~/Code/lib
 INCLUDEDIRS =  -I ~/Code/Integrator -I ~/Code/Tools  -I ~/Code/Glauber -I ~/Code/Knockout
STATICLIBS =  ~/Code/lib/libKnockout.a  ~/Code/lib/libGlauber.a  ~/Code/lib/libMePhys.a   ~/Code/lib/libIntegrator.a   


#		Header files
HDRS =		

# root headers

#		Object files
# OBJS = 	        $(EFILE:../../bin/=)  
# OBJS2 = 	        $(HDRS:.hpp=.o)  $(CPATH)/neutestint.o

#-------------------------------------------------------------------------
#		Executable file
EFILE   = 	../../bin/neutest ../../bin/neutestint
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
	$(CC) $(CFLAGS) $(@F).cpp $(LIBS) $(WIMLIBS) $(LDFLAGS) $(ROOT) $(LIBDIRS) $(INCLUDEDIRS) -o $@

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

