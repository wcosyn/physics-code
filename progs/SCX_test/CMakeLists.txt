project(SCX_test)

add_executable(SCX_test SCX_test.cpp )
target_link_libraries(SCX_test Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(SCX_test PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS SCX_test DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()