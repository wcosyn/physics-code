project(cuba_test)

add_executable(cubatest cubatest.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(cubatest Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(cubatest_wrapped cubatest_wrapped.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(cubatest_wrapped Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(cubatest_complex cubatest_complex.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(cubatest_complex Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(cubatest_multicomp cubatest_multicomp.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(cubatest_multicomp Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

if(NOT BUILD_PROGS)
set_target_properties(cubatest PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(cubatest_wrapped PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(cubatest_multicomp PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(cubatest_complex PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS cubatest cubatest_wrapped cubatest_multicomp cubatest_complex DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()