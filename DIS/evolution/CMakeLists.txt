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
add_library(Evolution SHARED ${EVOLUTION_SRCS})
target_include_directories(Evolution PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_link_libraries(Evolution Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)
set_target_properties(Evolution PROPERTIES ${EVOLUTION_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS Evolution DESTINATION "${CMAKE_SOURCE_DIR}/lib")

############### build the library #####################
if(STATIC_BUILD)
  add_library(EvolutionST STATIC ${EVOLUTION_SRCS})
  target_include_directories(EvolutionST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
 # target_link_libraries(EvolutionST Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a -Wl,-no-whole-archive)
  set_target_properties(EvolutionST PROPERTIES OUTPUT_NAME EVOLUTION)
  install(TARGETS EvolutionST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)