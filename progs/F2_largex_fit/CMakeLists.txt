project(F2largex)


#program used to fit the FSI parameters to the deeps data
add_executable(F2_fit F2_fit.cpp )  
target_link_libraries(F2_fit DIS Evolution MePhys)
add_executable(F2_bootstrap F2_bootstrap.cpp )  
target_link_libraries(F2_bootstrap DIS Evolution MePhys)

add_executable(EMC_fit EMC_fit.cpp)
target_link_libraries(EMC_fit MePhys)


if(NOT BUILD_PROGS)
set_target_properties(EMC_fit F2_fit F2_bootstrap PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS EMC_fit F2_fit F2_bootstrap DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()