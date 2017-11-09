project(pion)

add_executable(pion pion.cpp )
target_link_libraries(pion Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(pion PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS pion DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()