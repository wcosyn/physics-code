project(gpd)

add_executable(gpd gpd.cpp ) #small program to test he3 functionality for the LDRD project, never properly used so far.
target_link_libraries(gpd DIS MePhys)

if(NOT BUILD_PROGS)
set_target_properties(gpd PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS gpd DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()