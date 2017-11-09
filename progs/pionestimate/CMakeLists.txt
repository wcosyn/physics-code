project(pionestimate)

add_executable(pionestimate pionestimate.cpp )
target_link_libraries(pionestimate Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(pionestimate PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS pionestimate DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()