#originaly written by Florian Uhlig
#example of CMake script for compiling and linking root classse.  
#To create an Xcode project:
#   cmake -G Xcode ./ 
#extracted from FairRoot project (fairroot.gsi.de) by M. Al-Turany (m.al-turany@gsi.de)

#cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_BUILD_TYPE=Debug ..


# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.0 FATAL_ERROR)
MESSAGE("c++ ${CMAKE_CXX_COMPILER}")
MESSAGE("c ${CMAKE_C_COMPILER}")
MESSAGE("Type ${CMAKE_BUILD_TYPE}")
# Set name of our project to "Event". Has to be done
# after check of cmake version since this is a new feature                  
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
# Load some basic macros which are needed later on
find_package(ROOT REQUIRED)
find_package(GSL REQUIRED)
#Check the compiler and set the compile and link flags
# set(CMAKE_BUILD_TYPE RelWithDebInfo)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


# Create a library called "libEvent" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.
set(MPI_DIR "/opt/openmpi/1.4_intel/include")
set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} ${GSL_INCLUDE_DIRS} ${MPI_DIR}
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR} ${GSL_LIBRARY_DIRS}
)
 
link_directories( ${LINK_DIRECTORIES})
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})


enable_language (Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


EXEC_PROGRAM("hostname" OUTPUT_VARIABLE hostname)

if (hostname MATCHES accident)
   MESSAGE("--- Configuring for accident")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

#       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
#       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3  -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3 ")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g  -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3  -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3  -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O3 ")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O3  -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3  -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3  -fno-inline -ftest-coverage -fprofile-arcs")

      
      
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
      MESSAGE("--- Found Intel compiler collection")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0  -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3  -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3  -fno-inline -ftest-coverage -fprofile-arcs")
    endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")

  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	MESSAGE("--- Found Intel fortran compiler collection")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -coarray=shared   -shared-intel  -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel")
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



endif (hostname MATCHES accident) 

if (hostname MATCHES talvin13 OR hostname MATCHES trillian)
   MESSAGE("--- Configuring for talvin13")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

#       SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
#       SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3  -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3 ")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-O0 -g  -pg -fPIC")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O3  -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O3 ")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O3  -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3  -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-O0 -g  -pg -fPIC")

      
      
    endif (CMAKE_COMPILER_IS_GNUCXX)
    if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
      MESSAGE("--- Found Intel compiler collection")
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O0  -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O3  -shared-intel -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3  -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g -O0 -fno-reorder-blocks -fno-schedule-insns -fno-inline -prof_gen")
    endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")


  if (Fortran_COMPILER_NAME STREQUAL "gfortran")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
  elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
	MESSAGE("--- Found Intel fortran compiler collection")
    # ifort (untested)
    set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3 -coarray=shared   -shared-intel  -fp-model source -unroll -pthread -ipo -xsse4.2 -parallel")
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



endif (hostname MATCHES talvin13 OR hostname MATCHES trillian) 


set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR})
FIND_PROGRAM(XIAR xiar)
IF(XIAR)
  SET(CMAKE_AR "${XIAR}")
ENDIF(XIAR)
MARK_AS_ADVANCED(XIAR)


add_subdirectory(Tools)
add_subdirectory(Glauber)
add_subdirectory(Knockout)
add_subdirectory(DIS)

message(${CMAKE_CXX_FLAGS_RELEASE})
message(${CMAKE_Fortran_FLAGS_RELEASE})
message(${CMAKE_CXX_FLAGS_DEBUG})
message(${CMAKE_Fortran_FLAGS_DEBUG})