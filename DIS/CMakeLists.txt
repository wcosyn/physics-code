project(DIS)



# fill list of header files 
set(DIS_HEADERS_HPP 
InclusiveCross.hpp
DeuteronMomDistr.hpp
DeuteronCross.hpp
)

foreach(it ${DIS_HEADERS_HPP})
  STRING(REGEX REPLACE "(^.*).hpp" "\\1.cpp" bla ${it})
  set(DIS_SRCS ${DIS_SRCS} ${bla})
endforeach()


set(INCLUDE_DIRECTORIES
${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} "${CMAKE_SOURCE_DIR}/Tools" "${CMAKE_SOURCE_DIR}/Glauber"
)
include_directories( ${INCLUDE_DIRECTORIES})




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
target_link_libraries(DIS ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
set_target_properties(DIS PROPERTIES ${DIS_LIBRARY_PROPERTIES})
target_link_libraries (DIS MePhys Glauber)
############### build the library #####################
# add_library(DISST STATIC ${DIS_SRCS})
# # target_link_libraries(DIS ${ROOT_LIBRARIES} ${GSL_LIBRARIES})
# set_target_properties(DISST PROPERTIES OUTPUT_NAME DIS)

############### install the library ###################
install(TARGETS DIS LIBRARY DESTINATION "${CMAKE_SOURCE_DIR}/lib")
