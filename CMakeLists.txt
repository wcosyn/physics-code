#originaly written by Florian Uhlig
#example of CMake script for compiling and linking root classse.  
#To create an Xcode project:
#   cmake -G Xcode ./ 
#extracted from FairRoot project (fairroot.gsi.de) by M. Al-Turany (m.al-turany@gsi.de)

#cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_BUILD_TYPE=Debug ..


# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
project(mycode CXX Fortran)
MESSAGE("c++ ${CMAKE_CXX_COMPILER}")
MESSAGE("c ${CMAKE_C_COMPILER}")
MESSAGE("Type ${CMAKE_BUILD_TYPE}")
# Set name of our project to "Event". Has to be done
# after check of cmake version since this is a new feature                  
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# Load some basic macros which are needed later on
find_package(ROOT REQUIRED)
find_package(GSL REQUIRED)
# find_package(MPI REQUIRED)
find_package(Fortran REQUIRED)
find_package(Doxygen)
#Check the compiler and set the compile and link flags
# set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

set(HPC_NAMES 
haunter
gastly
dugtrio 
raichu
delcatty
gulpin
gligar
golett
phanpy
)


# Create a library called "libEvent" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.
set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} ${MPI_INCLUDE_PATH}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR} ${GSL_LIBRARY_DIRS}
)
 
link_directories( ${LINK_DIRECTORIES})
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})

set(LIBRARIES ${ROOT_LIBRARIES} "-lMathMore" ${GSL_LIBRARIES} ${MPI_LIBRARIES} ${Fortran_LIBRARIES})
# MESSAGE(${LIBRARIES})
enable_language (Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


EXEC_PROGRAM("hostname" OUTPUT_VARIABLE hostname)

if (hostname MATCHES accident OR hostname MATCHES talvin13 OR hostname MATCHES we05c807 OR hostname MATCHES we05c708)
   MESSAGE("--- Configuring for accident")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

#       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
#       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3  -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -mtune=native -std=c++11")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-O3 -pg ")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3  -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O3  -march=native -mtune=native")
      set(CMAKE_C_FLAGS_DEBUG            "-g")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3  -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-O3  -pg ")

      
      
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
      MESSAGE("--- Found Intel compiler collection")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0  -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3  -std=c++11 -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel -I/usr/include/x86_64-linux-gnu/c++/4.8")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-p -O0  -no-ipo")
    endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")

  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	MESSAGE("--- Found Intel fortran compiler collection")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -shared-intel  -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "g77")
    # g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  else (Fortran_COMPILER_NAME STREQUAL "gfortran")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif (Fortran_COMPILER_NAME STREQUAL "gfortran")



# endif (hostname MATCHES accident OR hostname MATCHES talvin13) 


# endif (hostname MATCHES trillian) 

#supercomputer config
elseif (hostname MATCHES .os)
foreach(hpcname ${HPC_NAMES})
if(hostname MATCHES ${hpcname})
   MESSAGE("--- Configuring for supercomputer")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

#       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
#       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3  -g -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -no-multibyte-chars")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3  -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -std=c++11 -no-multibyte-chars")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-O0 -g  -pg")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2  -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2 -no-multibyte-chars")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O3  -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3  -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-O0 -g  -pg")

      
      
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
      MESSAGE("--- Found Intel compiler collection")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -std=c++11 -shared-intel -fp-model source -unroll -pthread -ipo -xSSE4.1 -axSSE4.1,SSE4.2 -parallel -g -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -no-multibyte-chars")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -std=c++11 -shared-intel -fp-model source -unroll -pthread -ipo -xSSE4.1 -axSSE4.1,SSE4.2 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -no-multibyte-chars")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g  -std=c++11")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-p -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline -prof_gen")
      set(CMAKE_C_FLAGS_RELEASE          "-O3 -no-multibyte-chars")

    endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")


  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	MESSAGE("--- Found Intel fortran compiler collection")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -03 -shared-intel -fp-model source -unroll -ipo -xSSE4.1 -axSSE4.1,SSE4.2 -parallel -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "g77")
    # g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  else (Fortran_COMPILER_NAME STREQUAL "gfortran")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif (Fortran_COMPILER_NAME STREQUAL "gfortran")
endif(hostname MATCHES ${hpcname})
endforeach(hpcname ${HPC_NAMES})

else(hostname MATCHES accident OR hostname MATCHES talvin13 OR hostname MATCHES we05c807  OR hostname MATCHES we05c708) 
   MESSAGE("--- default configuration")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

#       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
#       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3  -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -mtune=native -std=c++11")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-O0 -pg ")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3  -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O3  -march=native -mtune=native")
      set(CMAKE_C_FLAGS_DEBUG            "-g")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3  -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-O0  -pg ")

      
      
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
      MESSAGE("--- Found Intel compiler collection")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0  -shared-intel -fp-model source -unroll -pthread -ipo -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3  -shared-intel -fp-model source -unroll -pthread -ipo")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-p -O0  -no-ipo")
    endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")

  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	MESSAGE("--- Found Intel fortran compiler collection")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3   -shared-intel  -fp-model source -unroll -pthread -ipo")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "g77")
    # g77
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
  else (Fortran_COMPILER_NAME STREQUAL "gfortran")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
  endif (Fortran_COMPILER_NAME STREQUAL "gfortran")


# endif (hostname MATCHES gengar OR hostname MATCHES haunter) 
endif (hostname MATCHES accident OR hostname MATCHES talvin13 OR hostname MATCHES we05c807  OR hostname MATCHES we05c708) 


set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR})
FIND_PROGRAM(XIAR xiar)
IF(XIAR)
  SET(CMAKE_AR "${XIAR}")
ENDIF(XIAR)
MARK_AS_ADVANCED(XIAR)


add_subdirectory(Integrator)
add_subdirectory(Tools)
add_subdirectory(Glauber)
add_subdirectory(Knockout)
add_subdirectory(DIS)

if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY IMMEDIATE)
  add_custom_target(doc
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  SOURCES ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif()

# message(${CMAKE_CXX_FLAGS_RELEASE})
# message(${CMAKE_Fortran_FLAGS_RELEASE})
# message(${CMAKE_CXX_FLAGS_DEBUG})
# message(${CMAKE_Fortran_FLAGS_DEBUG})
