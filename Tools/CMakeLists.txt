project(MePhys)



# fill list of header files 
set(MePhys_HEADERS_HPP 
NuclStructure.hpp
DeuteronStructure.hpp
Utilfunctions.hpp
LightConeKin2to2.hpp
LightConeKin2to3.hpp
He3wf.hpp
)

foreach(it ${MePhys_HEADERS_HPP})
  STRING(REGEX REPLACE "(^.*).hpp" "\\1.cpp" bla ${it})
  set(MePhys_SRCS ${MePhys_SRCS} ${bla})
endforeach()

set(MePhys_HEADERS_H 
TCstPWF.h
TDeuteron.h
TGrossPWF.h
TInterpolatingWavefunction.h
TWavefunctionImplementation.h
TYukawaPWF.h
Matrix.h
GammaStructure.h
FourVector.h
TSpinor.h
TLorentzQuaternion.h
TKinematics.h
TKinematics2to2.h
TElectronKinematics.h
pair.h
recmosh.h
newcoef.h
# TMPI.h
TKinematics2to3.h
TKinematics2to3WithLabAngles.h
TLeptonKinematics.h
mstwpdf.h
)
foreach(it ${MePhys_HEADERS_H})
  STRING(REGEX REPLACE "(^.*).h" "\\1.cpp" bla ${it})
  set(MePhys_SRCS ${MePhys_SRCS} ${bla})
endforeach()

set(MePhys_HEADERS ${MePhys_HEADERS_H} ${MePhys_HEADERS_HPP} 
# constants.hpp 
fortran.h
f77matrx.h
f77complx.h
f77char.h
)

# set everything needed for the root dictonary and create the
# dictionary
set(MePhys_LINKDEF 
TCstPWFLinkDef.h
TDeuteronLinkDef.h
TGrossPWFLinkDef.h
TInterpolatingWavefunctionLinkDef.h
TWavefunctionImplementationLinkDef.h
TYukawaPWFLinkDef.h
TSpinorLinkDef.h
TLorentzQuaternionLinkDef.h
TElectronKinematicsLinkDef.h
TKinematicsLinkDef.h
TKinematics2to2LinkDef.h
# TMPILinkDef.h
TKinematics2to3LinkDef.h
TKinematics2to3WithLabAnglesLinkDef.h
TLeptonKinematicsLinkDef.h
)


foreach(it ${MePhys_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dict.cpp" bla ${it})
  set(MePhys_DICTIONARY ${MePhys_DICTIONARY} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()

foreach(it ${MePhys_LINKDEF})
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dict_rdict.pcm" bla ${it})
  set(MePhys_PCM ${MePhys_PCM} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
endforeach()

# foreach(it ${MePhys_LINKDEF})
#   STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1Dictt.cpp" bla ${it})
#   set(MePhys_DICTIONARY_ST ${MePhys_DICTIONARY_ST} "${CMAKE_CURRENT_BINARY_DIR}/${bla}")
# endforeach()


set(MePhys_FORT 
alekhin.f
F1F209_test.f
cteq.f
spline.f
clebsch.f
getwavemod_tm.f
grsv2000pdf_g1.f
R199x.f
hmrs-b.f
)
# MESSAGE("${CMAKE_MODULE_PATH}")
# foreach(outfile ${MePhys_DICTIONARY}
# set(CURRENT_INCLUDE_DIRECTORIES
# ${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}
# )
# 
# include_directories( ${CURRENT_INCLUDE_DIRECTORIES})
set(INCLUDE_DIRECTORIES
${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
)
include_directories( ${INCLUDE_DIRECTORIES})

set(SHARE_DIR "${CMAKE_SOURCE_DIR}/share")

# message(${SHARE_DIR})

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp"
  )



list(LENGTH MePhys_DICTIONARY len1)
math(EXPR len2 "${len1} - 1")
foreach(val RANGE ${len2})
  list(GET MePhys_DICTIONARY ${val} val1)
#   list(GET MePhys_DICTIONARY_ST ${val} val3)
  list(GET MePhys_LINKDEF ${val} val2)
  STRING(REGEX REPLACE "(^.*)LinkDef.h" "\\1.h" bla ${val2})
  ROOT_GENERATE_DICTIONARY("${bla}" "${val2}" "${val1}" "${INCLUDE_DIRECTORIES}")
#   ROOT_GENERATE_DICTIONARY("${bla}" "${val2}" "${val3}" "${INCLUDE_DIRECTORIES}")
endforeach()



# add the dictionary to the list of source files
SET(MePhys_SRCS_LI ${MePhys_SRCS} ${MePhys_DICTIONARY} ${MePhys_FORT}) 
# SET(MePhys_SRCS_ST ${MePhys_SRCS} ${MePhys_DICTIONARY_ST} ${MePhys_FORT}) 


# Set the library version in the main CMakeLists.txt
SET(MePhys_MAJOR_VERSION 1)
SET(MePhys_MINOR_VERSION 0)
SET(MePhys_PATCH_VERSION 0)
SET(MePhys_VERSION "${MePhys_MAJOR_VERSION}.${MePhys_MINOR_VERSION}.${MePhys_PATCH_VERSION}")
SET(MePhys_LIBRARY_PROPERTIES ${MePhys_LIBRARY_PROPERTIES}
    VERSION "${MePhys_VERSION}"
    SOVERSION "${MePhys_MAJOR_VERSION}"
    SUFFIX ".so"
    CLEAN_DIRECT_OUTPUT 1
)


############### build the library #####################
add_library(MePhys SHARED ${MePhys_SRCS_LI})
target_link_libraries(MePhys ${LIBRARIES})
set_target_properties(MePhys PROPERTIES ${MePhys_LIBRARY_PROPERTIES})

############### install the library ###################
install(TARGETS MePhys DESTINATION "${CMAKE_SOURCE_DIR}/lib")
IF (60000 LESS found_vers)
foreach(it ${MePhys_PCM})
  install(FILES ${it}  DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endforeach()
ENDIF(60000 LESS found_vers)

############# build the library #####################
if(STATIC_BUILD)
  add_library(MePhysST STATIC ${MePhys_SRCS_LI})
  target_link_libraries(MePhysST -Wl,-no-whole-archive)
  set_target_properties(MePhysST PROPERTIES OUTPUT_NAME MePhys)
  install(TARGETS MePhysST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)

