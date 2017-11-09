project(inttest)

add_executable(inttest inttest.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(inttest Integrator)

if(NOT BUILD_PROGS)
set_target_properties(inttest PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS inttest DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()