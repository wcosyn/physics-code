project(DIS)



# fill list of header files 
set(DIS_SRCS 
InclusiveCross.cpp
# InclusiveCrossRes.cpp #excluded because it only used romberg integration routines and I don't use them any more (removed it from the library)
DeuteronCross.cpp
Resonance.cpp
He3Cross.cpp
Poldeut.cpp
GPD.cpp
TransGPD_set.cpp
)




# Set the library version in the main CMakeLists.txt
SET(DIS_MAJOR_VERSION 1)
SET(DIS_MINOR_VERSION 0)
SET(DIS_PATCH_VERSION 0)
SET(DIS_VERSION "${DIS_MAJOR_VERSION}.${DIS_MINOR_VERSION}.${DIS_PATCH_VERSION}")
SET(DIS_LIBRARY_PROPERTIES ${DIS_LIBRARY_PROPERTIES}
    VERSION "${DIS_VERSION}"
    SOVERSION "${DIS_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(DIS SHARED ${DIS_SRCS})
target_include_directories(DIS PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(DIS Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)
set_target_properties(DIS PROPERTIES ${DIS_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS DIS DESTINATION "${CMAKE_SOURCE_DIR}/lib")

############### build the library #####################
if(STATIC_BUILD)
  add_library(DISST STATIC ${DIS_SRCS})
  target_include_directories(DISST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
  target_link_libraries(DISST Knockout MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a -Wl,-no-whole-archive)
  set_target_properties(DISST PROPERTIES OUTPUT_NAME DIS)
  install(TARGETS DISST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)

add_subdirectory(evolution)