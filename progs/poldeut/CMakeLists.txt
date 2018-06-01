project(poldeut)

add_executable(poldeut poldeut.cpp )
target_link_libraries(poldeut DIS MePhys)

add_executable(wftest wftest.cpp )
target_link_libraries(wftest MePhys Integrator)

if(NOT BUILD_PROGS)
set_target_properties(poldeut wftest PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS poldeut wftest DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()