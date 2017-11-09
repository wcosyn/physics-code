project(inclusivedis)

add_executable(disincres disincres.cpp ) #program that computes F2D with FSI (3 resonance regions)
target_link_libraries(disincres DIS MePhys)

add_executable(f2n f2n.cpp ) #program that computes F2D with FSI (3 resonance regions)
target_link_libraries(f2n DIS MePhys)

add_executable(disuniform disuniform.cpp ) #program that computes F2D with FSI (3 resonance regions)
target_link_libraries(disuniform DIS MePhys)

add_executable(disincGauss disincGauss.cpp ) #program that computes F2D with FSI (3 resonance regions)
target_link_libraries(disincGauss DIS MePhys)

add_executable(Azzres Azzres.cpp ) #program that computes F2D with FSI (3 resonance regions)
target_link_libraries(Azzres DIS MePhys)

add_executable(Azzuniform Azzuniform.cpp ) #program that computes F2D with FSI (3 uniformonance regions)
target_link_libraries(Azzuniform DIS MePhys)

if(NOT BUILD_PROGS)
set_target_properties(disincres PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(f2n PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(disuniform PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(disincGauss PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(Azzres PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(Azzuniform PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS disincres f2n disuniform disincGauss Azzres Azzuniform DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()