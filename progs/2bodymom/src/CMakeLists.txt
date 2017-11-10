project(2body)

set(2b_SRCS 
    sampler_HALLA.cpp sampler_PRL99_072501.cpp sampler_HALLB.cpp parser.cpp event.cpp GlauberGridThick_SCX.cpp GlauberGridThick_SEL.cpp ClassGridThick_SCX.cpp)


#wrapper lib used in the python file MeanFieldNucleus.py in the same folder
add_library(MeanFieldNucleusWrapper SHARED MeanFieldNucleusWrapper.cpp)
target_include_directories(MeanFieldNucleusWrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(MeanFieldNucleusWrapper Glauber)
# set_target_properties(MeanFieldNucleusWrapper PROPERTIES ${MeanFieldNucleusWrapper_LIBRARY_PROPERTIES})

    
    
# all programs that camille used
add_executable(sampler sampler.cpp ${2b_SRCS})
target_link_libraries(sampler Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(2bodymom_new 2bodymom.cpp ${2b_SRCS})
target_link_libraries(2bodymom_new Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(SCX_test_new SCX_test.cpp ${2b_SRCS})
target_link_libraries(SCX_test_new Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(class_SCX class_SCX.cpp ${2b_SRCS})
target_link_libraries(class_SCX Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)

add_executable(sampler_test sampler_test.cpp ${2b_SRCS})
target_link_libraries(sampler_test Knockout Glauber MePhys ${CMAKE_SOURCE_DIR}/lib/libcuba.a)


if(NOT BUILD_PROGS)
set_target_properties(MeanFieldNucleusWrapper PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(sampler PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(sampler_test PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(SCX_test_new PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(class_SCX PROPERTIES EXCLUDE_FROM_ALL 1)
set_target_properties(2bodymom_new PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
install(TARGETS MeanFieldNucleusWrapper DESTINATION "${CMAKE_SOURCE_DIR}/lib")
install(TARGETS MeanFieldNucleusWrapper DESTINATION  ${CMAKE_CURRENT_SOURCE_DIR})
install(TARGETS sampler sampler_test SCX_test_new class_SCX 2bodymom_new DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()