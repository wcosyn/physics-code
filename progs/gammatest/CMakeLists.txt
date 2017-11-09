project(gamma)

add_executable(gamma gamma.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(gamma Knockout MePhys)

if(NOT BUILD_PROGS)
set_target_properties(gamma PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS gamma DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()