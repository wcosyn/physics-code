project(QE)

add_executable(transp transp.cpp )  #used to compute transparencies for Or for some proposal
target_link_libraries(transp Knockout Glauber MePhys)

add_executable(transp_jlab12 transp_jlab12.cpp )  #used to compute transparencies for Or for some proposal
target_link_libraries(transp_jlab12 Knockout Glauber MePhys)

add_executable(transp_uchannel transp_uchannel.cpp )  #used to compute transparencies for Or for some proposal
target_link_libraries(transp_uchannel Knockout Glauber MePhys)


add_executable(transp_eli transp_eli.cpp ) # used to compute lead transparencies for some proposal
target_link_libraries(transp_eli Knockout Glauber MePhys)

add_executable(A_QE_cross A_QE_cross.cpp )  #used to compute A(e,e'p) cross sections (Monoghan's paper and so)
target_link_libraries(A_QE_cross Knockout Glauber MePhys)

add_executable(JlabPAC16_calc JlabPAC16_calc.cpp )  #computed transparencies for the CAFE experiment proposal
target_link_libraries(JlabPAC16_calc Knockout Glauber MePhys)


if(NOT BUILD_PROGS)
set_target_properties(transp PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(transp_jlab12 PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(transp_eli PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(A_QE_cross PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(JlabPAC16_calc PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS transp_jlab12 transp transp_eli A_QE_cross JlabPAC16_calc DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()