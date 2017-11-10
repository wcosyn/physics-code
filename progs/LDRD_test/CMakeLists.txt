project(LDRD_test)

add_executable(LDRD_test LDRD_test.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(LDRD_test DIS MePhys)

if(NOT BUILD_PROGS)
set_target_properties(LDRD_test PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS LDRD_test DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()