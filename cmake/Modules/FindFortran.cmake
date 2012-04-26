# - Fortran System Library Module
#
# When linking C or C++ code to a library written in Fortran, it is necessary to
# specify the system's Fortran library.  This module attempts to detect the
# Fortran library in a system independent fashion. 
#
# This module will set the following variables:
#   Fortran_FOUND               TRUE if we found Fortran
#   Fortran_LIBRARIES           Link libraries for Fortran

include( FindPackageHandleStandardArgs )

# ensure that Fortran is enabled
get_property( _LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES )
if( NOT _LANGUAGES_ MATCHES Fortran )
  if( Fortran_FIND_REQUIRED )
    message( FATAL_ERROR 
      "Fortran must be enabled to interface C and Fortran code." )
  else()
    message( STATUS "Looking for Fortran... - NOT found (Fortran not enabled)")
    return()
  endif()
endif()

# The ld exe on GNU compilers is actually collect2
if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
    set( Fortran_LINKER_EXE collect2 )
else()
    get_filename_component( Fortran_LINKER_EXE ${CMAKE_LINKER} NAME )
endif()

message( STATUS "Fortran linker is ${Fortran_LINKER_EXE}" )

# compile a very simple test executable verbosely so the link command is
# output
get_filename_component( find_module_path ${CMAKE_CURRENT_LIST_FILE} PATH )

file( WRITE ${CMAKE_BINARY_DIR}/TestFortranCompileFlags.f
  "      program test\n"
  "      end\n"
)

set( compile_cmd
  ${CMAKE_Fortran_COMPILER} -v 
  ${CMAKE_BINARY_DIR}/TestFortranCompileFlags.f 
  -o TestFortranCompileFlags
)
execute_process( 
  COMMAND ${compile_cmd} 
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  RESULT_VARIABLE compile_result
  OUTPUT_VARIABLE compile_output
  ERROR_VARIABLE compile_output
)

if( NOT ${compile_result} EQUAL 0 )
  message( SEND_ERROR 
    "Unable to compile a simple Fortran test program.\n"
    "The return code or error message was:\n"
    "${compile_result}\n"
    "The compile output was: \n"
    "${compile_output}"
  )
endif()

# Extract any link commands from the compile output
string( REGEX MATCHALL "^.*${Fortran_LINKER_EXE}.*\$"
  link_cmd_list "${compile_output}"
)

# for each of the linker commands, strip out the -L flags to find the library
# paths.
foreach( link_cmd ${link_cmd_list} )
  # find all the -L... flags in the link line
  string( REGEX MATCHALL "-L([^\", ]+|\"[^\"]+\")" 
    Fortran_ALL_LINK_PATHS "${link_cmd}"
  )

  # strip the -L from the paths and add to the list of search locations
  set( Fortran_LIBRARY_SEARCH_PATH )
  foreach( lpath ${Fortran_ALL_LINK_PATHS} )
    string( REGEX REPLACE "^-L" "" lpath ${lpath} )
    string( REGEX REPLACE "//" "/" lpath ${lpath} )
    message( STATUS "Adding ${lpath} to Fortran library search path" )
    list( APPEND Fortran_LIBRARY_SEARCH_PATH ${lpath} )
  endforeach( lpath )
  
  # find all the -l... flags in the link line
  # match only -l's preceded by a space or comma
  # this is to exclude things like /redhat-linux
  string( REGEX MATCHALL "[, ]-l([^\", ]+)" 
    Fortran_ALL_LIBRARY_FLAGS "${link_cmd}"
  )
  
  # strip the leading -l from the link flags
  set( Fortran_ALL_LIBRARIES )
  foreach( lib ${Fortran_ALL_LIBRARY_FLAGS} )
    string( REGEX REPLACE "^[, ]-l" "" lib ${lib} )
    message( STATUS "Adding ${lib} to the list of Fortran libraries" )
    list( APPEND Fortran_ALL_LIBRARIES ${lib} )
  endforeach()
endforeach( link_cmd )

foreach( library ${Fortran_ALL_LIBRARIES} )
  find_library( Fortran_${library}_PATH 
    NAMES ${library} 
    PATHS ${Fortran_LIBRARY_SEARCH_PATH}
  )

  if( Fortran_${library}_PATH )
    list( APPEND Fortran_LIBRARIES ${Fortran_${library}_PATH} )
  endif()
endforeach()

find_package_handle_standard_args( Fortran DEFAULT_MSG Fortran_LIBRARIES )

