project(Evolution)



# fill list of header files 
set(EVOLUTION_SRCS 
decouple.f90
evolvers.f90
kinds.f90
msbar.f90
utils.f90
vector_nonsinglet.f90
vector_singlet.f90
wrappers.f90
)




# Set the library version in the main CMakeLists.txt
SET(EVOLUTION_MAJOR_VERSION 1)
SET(EVOLUTION_MINOR_VERSION 0)
SET(EVOLUTION_PATCH_VERSION 0)
SET(EVOLUTION_VERSION "${EVOLUTION_MAJOR_VERSION}.${EVOLUTION_MINOR_VERSION}.${EVOLUTION_PATCH_VERSION}")
SET(EVOLUTION_LIBRARY_PROPERTIES ${EVOLUTION_LIBRARY_PROPERTIES}
    VERSION "${EVOLUTION_VERSION}"
    SOVERSION "${EVOLUTION_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(EVOLUTION SHARED ${EVOLUTION_SRCS})
target_include_directories(EVOLUTION PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_link_libraries(EVOLUTION Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)
set_target_properties(EVOLUTION PROPERTIES ${EVOLUTION_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS EVOLUTION DESTINATION "${CMAKE_SOURCE_DIR}/lib")

############### build the library #####################
if(STATIC_BUILD)
  add_library(EVOLUTIONST STATIC ${EVOLUTION_SRCS})
  target_include_directories(EVOLUTIONST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
 # target_link_libraries(EVOLUTIONST Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a -Wl,-no-whole-archive)
  set_target_properties(EVOLUTIONST PROPERTIES OUTPUT_NAME EVOLUTION)
  install(TARGETS EVOLUTIONST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)