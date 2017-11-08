project(Glauber)



# fill list of source files 
set(Glauber_SRCS 
MeanFieldNucleus.cpp
MeanFieldNucleusThick.cpp
NucleonEMOperator.cpp
NucleonWeakOperator.cpp
FsiCorrelator.cpp
FastParticle.cpp
AbstractFsiGrid.cpp
AbstractFsiCTGrid.cpp
AbstractFsiCTDecayGrid.cpp
AbstractFsiGridThick.cpp
AbstractFsiDecayGridThick.cpp
AbstractFsiCTGridThick.cpp
AbstractFsiCTDecayGridThick.cpp
GlauberGrid.cpp
GlauberGridThick.cpp
GlauberDecayGridThick.cpp
OneGlauberGrid.cpp
TMFSpinor.cpp
DistMomDistrGrid.cpp
ROMEAGrid.cpp
ClassGridThick_SCX.cpp
)


# set everything needed for the root dictonary and create the
# dictionary
set(Glauber_LINKDEF 
TMFSpinorLinkDef.h
)

foreach(it ${Glauber_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dict_rdict.pcm" bla ${it})
  set(Glauber_PCM ${Glauber_PCM} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()


set(Glauber_FORT 
global2009.f
)


foreach(it ${Glauber_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dict.cpp" bla ${it})
  set(Glauber_DICTIONARY ${Glauber_DICTIONARY} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()

# set(INCLUDE_DIRECTORIES
# ${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_SOURCE_DIR}/Tools" "${CMAKE_SOURCE_DIR}/Integrator" "${CMAKE_BINARY_DIR}/Tools"
# )
# include_directories( ${INCLUDE_DIRECTORIES})

#this is a bit ugly but I don't see another way right now to make clear to the root generate dictionary function it depends on the MePhys lib too...
#Should probably throw out all these dict shit since noone will ever use it.
set(DicInc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/Tools)

list(LENGTH Glauber_DICTIONARY len1)
math(EXPR len2 "${len1} - 1")
foreach(val RANGE ${len2})
  list(GET Glauber_DICTIONARY ${val} val1)
  list(GET Glauber_LINKDEF ${val} val2)
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1.hpp" bla ${val2})
  ROOT_GENERATE_DICTIONARY("${bla}" "${val2}" "${val1}" "${DicInc}")
endforeach()

# add the dictionary to the list of source files
SET(Glauber_SRCS_LI ${Glauber_SRCS} ${Glauber_DICTIONARY} ${Glauber_FORT}) 


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
target_include_directories(Glauber PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Glauber MePhys Integrator)
set_target_properties(Glauber PROPERTIES ${Glauber_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS Glauber DESTINATION "${CMAKE_SOURCE_DIR}/lib")
IF(60000 LESS found_vers)
foreach(it ${Glauber_PCM})
  install(FILES ${it}  DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endforeach()
ENDIF(60000 LESS found_vers)
 
############## build the library #####################
if(STATIC_BUILD)
  add_library(GlauberST STATIC ${Glauber_SRCS_LI})
  target_include_directories(GlauberST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(GlauberST MePhys Integrator -Wl,-no-whole-archive)
  set_target_properties(GlauberST PROPERTIES OUTPUT_NAME Glauber)
  install(TARGETS GlauberST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)
