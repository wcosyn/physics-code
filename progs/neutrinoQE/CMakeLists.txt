project(neutrinoQE)

add_executable(ddiff ddiff.cpp )
target_link_libraries(ddiff Knockout Glauber MePhys)

add_executable(neutestfixEin neutestfixEin.cpp )
target_link_libraries(neutestfixEin Knockout Glauber MePhys)

add_executable(neutestint neutestint.cpp )
target_link_libraries(neutestint Knockout Glauber MePhys)



add_executable(sdiff_Q2 sdiff_Q2.cpp )
target_link_libraries(sdiff_Q2 Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(sdiff_Q2QE sdiff_Q2QE.cpp )
target_link_libraries(sdiff_Q2QE Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(sdiff_Q2p sdiff_Q2p.cpp )
target_link_libraries(sdiff_Q2p Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

if(NOT BUILD_PROGS)
set_target_properties(ddiff PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(neutestfixEin PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(neutestint PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(sdiff_Q2 PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(sdiff_Q2QE PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(sdiff_Q2p PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS ddiff neutestfixEin neutestint sdiff_Q2 sdiff_Q2QE sdiff_Q2p DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()