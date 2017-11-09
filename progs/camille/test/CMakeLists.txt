
include_directories(${PROJECT_SOURCE_DIR}/Integrator)
add_executable(test_numint integration_tests.cpp)
target_link_libraries(test_numint Integrator)

include_directories(${PROJECT_SOURCE_DIR}/Tools ${PROJECT_SOURCE_DIR}/Glauber)
add_executable(test_spinor spinor_test.cpp)
target_link_libraries(test_spinor Glauber MePhys Integrator)

add_executable(test_glauber glauber_test.cpp)
target_link_libraries(test_glauber Glauber MePhys Integrator)

include_directories(${PROJECT_SOURCE_DIR}/Knockout)
add_executable(test_matrixEl matrixEl_test.cpp)
target_link_libraries(test_matrixEl MePhys Glauber Integrator Knockout)
set_target_properties(test_matrixEl PROPERTIES COMPILE_FLAGS "-Wall")


add_executable(test_tmfspinor tmfspinor_test.cpp)
target_link_libraries(test_tmfspinor MePhys Glauber Integrator Knockout)

add_executable(test_tmfspinor2D tmfspinor2D_test.cpp)
target_link_libraries(test_tmfspinor2D MePhys Glauber Integrator Knockout)

if(NOT BUILD_PROGS)
set_target_properties(test_spinor PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(test_glauber PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(test_matrixEl PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(test_tmfspinor PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(test_tmfspinor2D PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install( TARGETS test_spinor test_glauber test_matrixEl test_tmfspinor test_tmfspinor2D
	RUNTIME DESTINATION  ${CMAKE_SOURCE_DIR}//bin)
endif()

# add_custom_command(TARGET test DEPENDS tmfspinor2D_test.cpp
# 	COMMAND ${CMAKE_COMMAND} add_executable(test_tmfspinor2D tmfspinor2D_test.cpp)
# 	COMMAND ${CMAKE_COMMAND} target_link_libraries(test_tmfspinor2D MePhys Glauber Integrator Knockout)
# 	COMMAND ${CMAKE_COMMAND} install(TARGETS test_tmfspinor2D RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin/test)
# 	)

# add_custom_target(test DEPENDS test_tmfspinor2D COMMAND ${CMAKE_COMMAND} "test")
