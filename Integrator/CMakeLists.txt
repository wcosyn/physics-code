project(INTEGRATOR)

set(INTEGRATOR_HEADERS
numint/array.hpp
numint/macros.hpp
numint/numint.hpp
numint/numint_tools.hpp
numint/quadrature.hpp
numint/typedef.hpp
numint/numint2adapt.hpp
adaptive/cubature.h)


set(INTEGRATOR_SRCS ${INTEGRATOR_SRCS}
numint/cube_romb.cpp
numint/cube_tr.cpp
numint/fubini.cpp
numint/quad_cc.cpp
numint/quad_romb.cpp
numint/quad_tr.cpp
numint/vfubini.cpp
adaptive/cubature.c
)

set(INCLUDE_DIRECTORIES
${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/numint )
include_directories( ${INCLUDE_DIRECTORIES})




# Set the library version in the main CMakeLists.txt
SET(INTEGRATOR_MAJOR_VERSION 1)
SET(INTEGRATOR_MINOR_VERSION 0)
SET(INTEGRATOR_PATCH_VERSION 0)
SET(INTEGRATOR_VERSION "${INTEGRATOR_MAJOR_VERSION}.${INTEGRATOR_MINOR_VERSION}.${INTEGRATOR_PATCH_VERSION}")
SET(INTEGRATOR_LIBRARY_PROPERTIES ${INTEGRATOR_LIBRARY_PROPERTIES}
    VERSION "${INTEGRATOR_VERSION}"
    SOVERSION "${INTEGRATOR_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(Integrator SHARED ${INTEGRATOR_SRCS})
target_link_libraries(Integrator ${LIBRARIES})
set_target_properties(Integrator PROPERTIES ${INTEGRATOR_LIBRARY_PROPERTIES})
# target_link_libraries (INTEGRATOR MePhys Glauber Knockout)
############### build the library #####################
add_library(IntegratorST STATIC ${INTEGRATOR_SRCS})
target_link_libraries(IntegratorST -Wl,-no-whole-archive)
set_target_properties(IntegratorST PROPERTIES OUTPUT_NAME Integrator)

############### install the library ###################
install(TARGETS Integrator DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS IntegratorST DESTINATION "${CMAKE_SOURCE_DIR}/lib")