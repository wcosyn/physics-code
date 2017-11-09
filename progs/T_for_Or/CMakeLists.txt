project(transp)

add_executable(transp transp.cpp )
target_link_libraries(transp Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(transp PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS transp DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()