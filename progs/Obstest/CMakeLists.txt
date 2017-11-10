project(obstest)

add_executable(observables observables.cpp )
target_link_libraries(observables Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(observables PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS observables DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()