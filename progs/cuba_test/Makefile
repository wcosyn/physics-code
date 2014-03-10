# Makefile with auto dep checking
# Note that the build will fail if header files are renamed
# make will crash reporting it can't find the renamed header files
# because the *.d depency files will be out of date.
# This can be fixed with extra code but to keep the Makefile nice & short
# Just do make clean and make to rebuild the correct dependency files
#
# @author Camille Colle


CC     = g++
EXTRAFLAGS = #-Wall -Wextra -Weffc++
CFLAGS = -O3 $(EXTRAFLAGS) `gsl-config --cflags` `root-config --cflags`
INC    = -I$(HOME)/Code/Tools -I$(HOME)/Code/Integrator
LIBS   = -L$(HOME)/Code/lib -lcuba -lIntegrator `gsl-config --libs` `root-config --libs`


TARGETS=cubatest cubatest_wrapped cubatest_multicomp # list here all files that contain a main()
OBJS = cubatest.o cubatest_wrapped.o cubatest_multicomp.o # list here all object files needed to build
DEPS:=$(OBJS:.o=.d)

all: $(TARGETS)

$(TARGETS): %: %.o
	@tput setf 4
	@tput bold
	@printf "Building executable $@ \n"
	@tput sgr0
	$(CC) $^ $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"
-include $(DEPS)


%.o: %.cpp
	@tput setf 6
	@printf "Building dependencies for $@ \n"
	@tput sgr0
	$(CC) $< -MM $(INC) > $(patsubst %.o,%.d,$@)
	@tput setf 6
	@printf "Building object file $@ \n"
	@tput sgr0
	$(CC) $< -c $(CFLAGS) $(INC) -o $@
	@printf "\n"
.PHONY: clean
clean:
	rm *.o *.d $(TARGETS)
