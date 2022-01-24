project(gpd_progs)

include_directories(${QT_INCLUDE_DIRS} ${SFML_INCLUDE_DIR} ${ElementaryUtils_INCLUDE_DIR} ${NumA++_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR})

add_executable(conv_T_test conv_T_test.cpp ) #tranversity GPD program
target_link_libraries(conv_T_test GPD MePhys)

add_executable(two_vector_N two_vector_N.cpp)
target_link_libraries(two_vector_N GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(diDVCS diDVCS.cpp)
target_link_libraries(diDVCS GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(diDVCS_Q2 diDVCS_Q2.cpp)
target_link_libraries(diDVCS_Q2 GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(diDVCS_gen diDVCS_gen.cpp)
target_link_libraries(diDVCS_gen GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(diDVCS_gen_grid diDVCS_gen_grid.cpp)
target_link_libraries(diDVCS_gen_grid GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(diDVCS_gen_gridonly diDVCS_gen_gridonly.cpp)
target_link_libraries(diDVCS_gen_gridonly ${ROOT_LIBRARIES})

add_executable(two_vector_electro_N two_vector_electro_N.cpp)
target_link_libraries(two_vector_electro_N GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(two_vector_D two_vector_D.cpp)
target_link_libraries(two_vector_D GPD   ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(conv_V_test conv_V_test.cpp)
target_link_libraries(conv_V_test GPD  MePhys Glauber ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})

add_executable(Nucl_GPD Nucl_GPD.cpp)
target_link_libraries(Nucl_GPD GPD ${ElementaryUtils_LIBRARIES} ${QT_LIBRARIES})


if(NOT BUILD_PROGS)
set_target_properties(conv_T_test PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(two_vector_N PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(two_vector_electro_N PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(two_vector_D PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(conv_V_test PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(Nucl_GPD PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(diDVCS PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(diDVCS_Q2 PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(diDVCS_gen PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(diDVCS_gen_grid PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(diDVCS_gen_gridonly PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS conv_T_test DESTINATION "${CMAKE_SOURCE_DIR}/bin")
install(TARGETS two_vector_N DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS two_vector_electro_N DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS two_vector_D DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS conv_V_test DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS Nucl_GPD DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS diDVCS DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS diDVCS_Q2 DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS diDVCS_gen DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS diDVCS_gen_grid DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
install(TARGETS diDVCS_gen_gridonly DESTINATION "${CMAKE_SOURCE_DIR}/GPD_Partons/bin")
endif()
