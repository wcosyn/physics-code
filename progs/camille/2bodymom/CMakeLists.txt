project(2bodymom_run)

add_executable(2bodymom 2bodymom.cpp ) #camille's prog, don't know
target_link_libraries(2bodymom Knockout)

add_executable(2bodymom_glauber 2bodymom_glauber.cpp ) #camille's prog, don't know
target_link_libraries(2bodymom_glauber Knockout)

if(NOT BUILD_PROGS)
set_target_properties(2bodymom PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(2bodymom_glauber PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS 2bodymom 2bodymom_glauber DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()