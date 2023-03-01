project(obstest)

add_executable(observables observables.cpp )
target_link_libraries(observables Knockout Glauber MePhys)
add_executable(Ca_dens Ca_dens.cpp )
target_link_libraries(Ca_dens Knockout Glauber MePhys)

add_executable(Nucl_dens_profile Nucl_dens_profile.cpp )
target_link_libraries(Nucl_dens_profile Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
set_target_properties(observables PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS observables Ca_dens Nucl_dens_profile DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()