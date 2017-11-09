project(DQEinc)

add_executable(DQEinc DQEinc.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(DQEinc Knockout)

if(NOT BUILD_PROGS)
set_target_properties(DQEinc PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS DQEinc DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()