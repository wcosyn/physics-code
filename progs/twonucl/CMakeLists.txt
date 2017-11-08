project(twonucl)

add_executable(twonucl twonucl.cpp )
target_link_libraries(twonucl Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(twonucl PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS twonucl DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()