project(RhoT)


if(MPI_ENABLED)

    add_executable(rhottest rhottest.cpp )  #program used to test rho transparency calculations, I think never used for proper calculations
    target_link_libraries(rhottest Knockout Glauber MePhys)

    add_executable(rhod rhod.cpp )  #deuteron rho cross section, never used for actual calculations
    target_link_libraries(rhod Knockout Glauber MePhys)

    if(NOT BUILD_PROGS)
        set_target_properties(rhottest PROPERTIES EXCLUDE_FROM_ALL 1)
        set_target_properties(rhod PROPERTIES EXCLUDE_FROM_ALL 1)
        endif()

    if(BUILD_PROGS)
        install(TARGETS rhottest rhod DESTINATION "${CMAKE_SOURCE_DIR}/bin")
    endif()
endif()

add_executable(rhotest rhotest.cpp )  #program used to test rho grids (glauber and momentum distributions)
target_link_libraries(rhotest Knockout Glauber MePhys)

add_executable(rhosingle rhosingle.cpp )  #program used to compute single node rho transparencies, never used for actual calculations
target_link_libraries(rhosingle Knockout Glauber MePhys)

add_executable(rhoset rhoset.cpp )  #program used to compute single node rho transparencies, with the sigma in glauber set.  This was the program used for the actual calculations
target_link_libraries(rhoset Knockout Glauber MePhys)

add_executable(rhodint rhodint.cpp )  #program used for single node deuteron rho cross sections, never used in actual calculations
target_link_libraries(rhodint Knockout Glauber MePhys)

if(NOT BUILD_PROGS)
    set_target_properties(rhotest PROPERTIES EXCLUDE_FROM_ALL 1)
    set_target_properties(rhosingle PROPERTIES EXCLUDE_FROM_ALL 1)
    set_target_properties(rhoset PROPERTIES EXCLUDE_FROM_ALL 1)
    set_target_properties(rhodint PROPERTIES EXCLUDE_FROM_ALL 1)
    endif()

if(BUILD_PROGS)
    install(TARGETS rhotest rhosingle rhoset rhodint DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()
