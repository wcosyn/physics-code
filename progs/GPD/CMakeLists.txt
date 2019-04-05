project(gpd_progs)

include_directories(${QT_INCLUDE_DIRS} ${SFML_INCLUDE_DIR} ${ELEMENTARY_UTILS_INCLUDE_DIR} ${NUMA_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR})

add_executable(conv_T_test conv_T_test.cpp ) #tranversity GPD program
target_link_libraries(conv_T_test GPD MePhys)

add_executable(two_vector_N two_vector_N.cpp)
target_link_libraries(two_vector_N GPD   ${ELEMENTARY_UTILS_LIBRARIES} ${QT_LIBRARIES})

add_executable(conv_V_test conv_V_test.cpp)
target_link_libraries(conv_V_test GPD  MePhys ${ELEMENTARY_UTILS_LIBRARIES} ${QT_LIBRARIES})

if(NOT BUILD_PROGS)
set_target_properties(conv_T_test PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(two_vector_N PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(conv_V_test PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS conv_T_test DESTINATION "${CMAKE_SOURCE_DIR}/bin")
install(TARGETS two_vector_N DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS conv_V_test DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
endif()
