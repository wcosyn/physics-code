project(gpd_progs)

include_directories(${QT_INCLUDE_DIRS} ${SFML_INCLUDE_DIR} ${ELEMENTARY_UTILS_INCLUDE_DIR} ${NUMA_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR})

add_executable(gpd gpd.cpp ) #tranversity GPD program
target_link_libraries(gpd GPD MePhys)

add_executable(two_vector_N two_vector_N.cpp)
target_link_libraries(two_vector_N GPD   ${ELEMENTARY_UTILS_LIBRARIES} ${QT_LIBRARIES})

if(NOT BUILD_PROGS)
set_target_properties(gpd PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(two_vector_N PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS gpd DESTINATION "${CMAKE_SOURCE_DIR}/bin")
install(TARGETS two_vector_N DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
endif()
