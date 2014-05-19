project(Glauber)



# fill list of header files 
set(Glauber_HEADERS_HPP 
MeanFieldNucleus.hpp
MeanFieldNucleusThick.hpp
NucleonEMOperator.hpp
NucleonWeakOperator.hpp
FsiCorrelator.hpp
FastParticle.hpp
AbstractFsiGrid.hpp
AbstractFsiCTGrid.hpp
AbstractFsiCTDecayGrid.hpp
AbstractFsiGridThick.hpp
AbstractFsiDecayGridThick.hpp
AbstractFsiCTGridThick.hpp
AbstractFsiCTDecayGridThick.hpp
GlauberGrid.hpp
GlauberGridThick.hpp
GlauberDecayGridThick.hpp
OneGlauberGrid.hpp
TMFSpinor.hpp
DistMomDistrGrid.hpp
ROMEAGrid.hpp
)

foreach(it ${Glauber_HEADERS_HPP})
  STRING(REGEX REPLACE "(^.*).hpp" "\\1.cpp" bla ${it})
  set(Glauber_SRCS ${Glauber_SRCS} ${bla})
endforeach()

set(Glauber_HEADERS_H 
)
foreach(it ${Glauber_HEADERS_H})
  STRING(REGEX REPLACE "(^.*).h" "\\1.cpp" bla ${it})
  set(Glauber_SRCS ${Glauber_SRCS} ${bla})
endforeach()

set(Glauber_HEADERS ${Glauber_HEADERS_H} ${Glauber_HEADERS_HPP})

# set everything needed for the root dictonary and create the
# dictionary
set(Glauber_LINKDEF 
TMFSpinorLinkDef.h
)

set(Glauber_FORT 
global2009.f
)


foreach(it ${Glauber_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dict.cpp" bla ${it})
  set(Glauber_DICTIONARY ${Glauber_DICTIONARY} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()

foreach(it ${Glauber_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dictt.cpp" bla ${it})
  set(Glauber_DICTIONARY_ST ${Glauber_DICTIONARY_ST} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()

set(INCLUDE_DIRECTORIES
${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_SOURCE_DIR}/Tools" "${CMAKE_SOURCE_DIR}/Integrator"
)
include_directories( ${INCLUDE_DIRECTORIES})

list(LENGTH Glauber_DICTIONARY len1)
math(EXPR len2 "${len1} - 1")
foreach(val RANGE ${len2})
  list(GET Glauber_DICTIONARY ${val} val1)
  list(GET Glauber_DICTIONARY_ST ${val} val3)
  list(GET Glauber_LINKDEF ${val} val2)
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1.hpp" bla ${val2})
  ROOT_GENERATE_DICTIONARY("${bla}" "${val2}" "${val1}" "${INCLUDE_DIRECTORIES}")
  ROOT_GENERATE_DICTIONARY("${bla}" "${val2}" "${val3}" "${INCLUDE_DIRECTORIES}")
endforeach()

# add the dictionary to the list of source files
SET(Glauber_SRCS_LI ${Glauber_SRCS} ${Glauber_DICTIONARY} ${Glauber_FORT}) 
SET(Glauber_SRCS_ST ${Glauber_SRCS} ${Glauber_DICTIONARY_ST} ${Glauber_FORT}) 


# Set the library version in the main CMakeLists.txt
SET(Glauber_MAJOR_VERSION 1)
SET(Glauber_MINOR_VERSION 0)
SET(Glauber_PATCH_VERSION 0)
SET(Glauber_VERSION "${Glauber_MAJOR_VERSION}.${Glauber_MINOR_VERSION}.${Glauber_PATCH_VERSION}")
SET(Glauber_LIBRARY_PROPERTIES ${Glauber_LIBRARY_PROPERTIES}
    VERSION "${Glauber_VERSION}"
    SOVERSION "${Glauber_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(Glauber SHARED ${Glauber_SRCS_LI})
target_link_libraries(Glauber ${LIBRARIES} MePhys Integrator)
 set_target_properties(Glauber PROPERTIES ${Glauber_LIBRARY_PROPERTIES})
############## build the library #####################
add_library(GlauberST STATIC ${Glauber_SRCS_ST})
target_link_libraries(GlauberST -Wl,-no-whole-archive)
set_target_properties(GlauberST PROPERTIES OUTPUT_NAME Glauber)

############### install the library ###################
install(TARGETS Glauber DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS GlauberST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
