project(disres)

add_executable(increstest increstest.cpp )
target_link_libraries(increstest DIS MePhys)

add_executable(threeres threeres.cpp )
target_link_libraries(threeres DIS MePhys)

if(NOT BUILD_PROGS)
set_target_properties(increstest PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(threeres PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS increstest threeres DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()