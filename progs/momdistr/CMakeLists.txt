project(momdistr)

add_executable(momdistr momdistr.cpp )
target_link_libraries(momdistr Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(momdistr PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS momdistr DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()