include_directories(${PROJECT_SOURCE_DIR}/Tools ${PROJECT_SOURCE_DIR}/Knockout ${PROJECT_SOURCE_DIR}/Glauber)
add_executable(2bodycom 2bodycom.cpp)
target_link_libraries(2bodycom MePhys Knockout Glauber Integrator)
install(TARGETS 2bodycom
	RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin)
