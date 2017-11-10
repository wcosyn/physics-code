project(Knockout)



# fill list of header files 
set(Knockout_SRCS 
DeuteronMomDistr.cpp
Model.cpp
Cross.cpp
DoubleNModel.cpp
DoubleNCross.cpp
RhoTCross.cpp
RhoDeuteron.cpp
WeakQEHadronCurrent.cpp
WeakQECross.cpp
DQEinclusive.cpp
)




# Set the library version in the main CMakeLists.txt
SET(Knockout_MAJOR_VERSION 1)
SET(Knockout_MINOR_VERSION 0)
SET(Knockout_PATCH_VERSION 0)
SET(Knockout_VERSION "${Knockout_MAJOR_VERSION}.${Knockout_MINOR_VERSION}.${Knockout_PATCH_VERSION}")
SET(Knockout_LIBRARY_PROPERTIES ${Knockout_LIBRARY_PROPERTIES}
    VERSION "${Knockout_VERSION}"
    SOVERSION "${Knockout_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(Knockout SHARED ${Knockout_SRCS})
target_include_directories(Knockout PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Knockout MePhys Integrator Glauber)
set_target_properties(Knockout PROPERTIES ${Knockout_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS Knockout DESTINATION "${CMAKE_SOURCE_DIR}/lib")


############### build the library #####################
if(STATIC_BUILD)
  add_library(KnockoutST STATIC ${Knockout_SRCS})
  target_include_directories(KnockoutST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(KnockoutST MePhys Integrator Glauber -Wl,-no-whole-archive)
  set_target_properties(KnockoutST PROPERTIES OUTPUT_NAME Knockout)
  install(TARGETS KnockoutST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)
