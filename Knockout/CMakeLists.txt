project(Knockout)



# fill list of header files 
set(Knockout_HEADERS_HPP 
DeuteronMomDistr.hpp
)

foreach(it ${Knockout_HEADERS_HPP})
  STRING(REGEX REPLACE "(^.*).hpp" "\\1.cpp" bla ${it})
  set(Knockout_SRCS ${Knockout_SRCS} ${bla})
endforeach()


set(INCLUDE_DIRECTORIES
${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_SOURCE_DIR}/Integrator" 
"${CMAKE_SOURCE_DIR}/Tools" "${CMAKE_SOURCE_DIR}/Glauber"
)
include_directories( ${INCLUDE_DIRECTORIES})




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
target_link_libraries(Knockout ${LIBRARIES} MePhys Integrator)
set_target_properties(Knockout PROPERTIES ${Knockout_LIBRARY_PROPERTIES})
############### build the library #####################
add_library(KnockoutST STATIC ${Knockout_SRCS})
target_link_libraries(KnockoutST -Wl,-no-whole-archive)
set_target_properties(KnockoutST PROPERTIES OUTPUT_NAME Knockout)

############### install the library ###################
install(TARGETS Knockout DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS KnockoutST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
