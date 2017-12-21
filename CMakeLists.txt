#originaly written by Florian Uhlig
#example of CMake script for compiling and linking root classse.  
#To create an Xcode project:
#   cmake -G Xcode ./ 
#extracted from FairRoot project (fairroot.gsi.de) by M. Al-Turany (m.al-turany@gsi.de)

#cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_BUILD_TYPE=Debug ..


# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set defaults
if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT BUILD_PROGS)
  set(BUILD_PROGS no)
endif()

if(NOT STATIC_BUILD)
  set(STATIC_BUILD no)
endif()

if(NOT MPI_ENABLED)
  set(MPI_ENABLED no)
endif()

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
if(MPI_ENABLED)
  find_package(MPI REQUIRED)
  MESSAGE("MPI compile flags ${MPI_COMPILE_FLAGS}" )
  MESSAGE("MPI link flags ${MPI_LINK_FLAGS}" )
  MESSAGE("MPI lib ${MPI_LIBRARY}" )
  MESSAGE("MPI lib alt ${MPI_EXTRA_LIBRARY}" )
  MESSAGE("MPI libs ${MPI_LIBRARIES}" )
  endif()

find_package(Fortran REQUIRED)
find_package(Doxygen)

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) #to create compile_commands.json file used in VS code

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


 
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR})
set(CMAKE_INSTALL_RPATH "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH 1)

enable_language (Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_COMPILER_IS_CLANGCXX 0)
set(CMAKE_COMPILER_IS_CLANGCC 0)
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")	 	
  set(CMAKE_COMPILER_IS_GNUCXX 1)
  set(CMAKE_COMPILER_IS_CLANGCXX 1)
endif()
if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_COMPILER_IS_GNUCC 1)
  set(CMAKE_COMPILER_IS_CLANGCC 1)
endif()

EXEC_PROGRAM("hostname" OUTPUT_VARIABLE hostname)

set(RELEASE_FLAGS -Wall -Wextra)
set(DEBUG_FLAGS ${RELEASE_FLAGS})

if (CMAKE_COMPILER_IS_GNUCXX)
MESSAGE("--- Found GNU compiler collection")   
  set(RELEASE_FLAGS ${RELEASE_FLAGS} -O3 -march=native -mtune=native -pipe)
  set(DEBUG_FLAGS ${DEBUG_FLAGS} -g)
endif (CMAKE_COMPILER_IS_GNUCXX)
  
if(${CMAKE_CXX_COMPILER} MATCHES "icpc")
  MESSAGE("--- Found Intel compiler collection")
  set(RELEASE_FLAGS ${RELEASE_FLAGS} -O3 -shared-intel -fp-model source -unroll -pthread -ipo -xavx)
  set(DEBUG_FLAGS ${DEBUG_FLAGS} -g)

  if(MPI_ENABLED)
   MESSAGE("-- Enabling mpi flags")
  set(RELEASE_FLAGS ${RELEASE_FLAGS} -parallel)  
  endif()  
endif(${CMAKE_CXX_COMPILER} MATCHES "icpc")


# if (hostname MATCHES we05c807 OR hostname MATCHES we05c708 OR hostname MATCHES we05c700)
#    MESSAGE("--- Configuring for machine in our group")

#   if (Fortran_COMPILER_NAME STREQUAL "gfortran")
#   elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
#   elseif (Fortran_COMPILER_NAME STREQUAL "g77")
#   else (Fortran_COMPILER_NAME STREQUAL "gfortran")
#     message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
#     message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
#     message ("No optimized Fortran compiler flags are known, we just try -O2...")
#   endif (Fortran_COMPILER_NAME STREQUAL "gfortran")
# #supercomputer config
# elseif (hostname MATCHES .os)
# foreach(hpcname ${HPC_NAMES})
# if(hostname MATCHES ${hpcname})
#    MESSAGE("--- Configuring for supercomputer")
# endif(hostname MATCHES ${hpcname})
# endforeach(hpcname ${HPC_NAMES})

# else(hostname MATCHES we05c807  OR hostname MATCHES we05c708 OR hostname MATCHES we05c700) 
#    MESSAGE("--- default configuration")
# endif (hostname MATCHES we05c807  OR hostname MATCHES we05c708 OR hostname MATCHES we05c700 ) 


# set(CMAKE_PREFIX_PATH ${CMAKE_SOURCE_DIR})
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
add_subdirectory(progs)

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
