# A Makefile with auto dep generation
# Note that the build will fail if header
# files are renamed. This is because the
# *.d dependency files will be out of date
# and point to non-existing header files
# We can solve this automatically but to
# keep the code nice and simple just do
# $> make clean (or rm *.d)
# $> make
# If header files have been renamed
#
# @author Camille Colle

CC=g++
WFLAGS=-Wall -Wextra -Wshadow #-Weffc++ # You can expect to see a buttload of warnings when enabling all these flags
CFLAGS= -O2 -std=c++11 -fPIC $(WFLAGS) `gsl-config --cflags` `root-config --cflags`
# -fext-numeric-literals is there to allow conversion between C and C++ complex numbers
#
# -Wl means pass arguments to runtime linker. -Wl is followed by comma seperated list with arguments to be passed to runtime linker
#  This is needed because the -L, -l flags only link @compile time but not runtime! If the libraries are not in the standard library paths
#  or LD_LIBRARY_PATH is not set we will get in trouble because the dynamic linker will (ld) will not be able to find our libs @runtime. 
#  To avoid using LD_LIBRARY_PATH we use rpath...
LIBS=-L$(HOME)/Code/trunk/lib -Wl,-rpath=$(HOME)/Code/trunk/lib -lcuba -lIntegrator -lGlauber -lKnockout -lMePhys `gsl-config --libs` `root-config --libs`
INC=-I$(HOME)/Code/trunk/Tools -I$(HOME)/Code/trunk/Glauber -I$(HOME)/Code/trunk/Integrator

TARGETS=sampler 2bodymom SCX_test RPWIA_toy#sampler_test       # files with a main()
OBJS=sampler_HALLA.o sampler_PRL99_072501.o sampler_HALLB.o parser.o event.o GlauberGridThick_SCX.o GlauberGridThick_SEL.o # object files for building, no main()'s here !!!
SOLIBS=MeanFieldNucleusWrapper.so
DEPDIR:=dep
DEPS=$(patsubst %.o, $(DEPDIR)/%.d, $(OBJS))

.DEFAULT: all
all: $(TARGETS) $(SOLIBS)


-include $(DEPS)


# if object_deps gets called from rule with target $@ will contain the target's name
define object_deps
	@tput setf 3; tput bold
	@printf "\n\nobject_deps called with arguments $@\n\n";
	@cat $(patsubs %,%.cpp,$@) | grep -E "\<.*\>"
	@tput sgr0;
endef

$(SOLIBS): %.so: %.o
	@tput setaf 4; tput bold
	@printf "Building shared object $@ : $^ \n"
	@tput sgr0
	$(CC) -shared $< $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"

# dependency of each target is all the object files and it's own %.o file.
# Technically too much dependencies... (targets will generally not depend on all $(OBJS) files)
# But it's better to have too much dependencies then too few.
$(TARGETS): %: %.o $(OBJS)
	@tput setaf 1; tput bold
	@printf "Linking executable $@  :  $^ \n"
	@tput sgr0
	$(CC) $^ $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"
	#mv $@ ../run/


%.o : %.cpp
	@echo DEPS INCLUDED ARE $(DEPS)
	@mkdir -p $(DEPDIR)
	@tput setf 6
	@printf "Building dependencies for $<\n"
	@tput sgr0
	$(CC) -MM $< $(CFLAGS) $(INC) > $(patsubst %.o,$(DEPDIR)/%.d,$@)
	@tput setf 6
	@printf "Building  $@ :  $^ \n"
	@tput sgr0
	$(CC) -c $< $(CFLAGS) $(INC) -o $@

.PHONY: clean

clean:
	rm *.o $(DEPDIR)/*.d *.so 
