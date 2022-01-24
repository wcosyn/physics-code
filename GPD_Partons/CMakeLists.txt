project(GPD)

# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
#cmake_minimum_required(VERSION 2.6)

# define project name and its language
#project(PARTONS_release CXX)

# define c++ standard and issue all the warning demanded by this standard
#add_definitions(-std=c++98 -pedantic)

# AUTOMOC ==================================================================================

# tell cmake to run Qt moc when necessary
#set(CMAKE_AUTOMOC ON)
  
# as moc files are generated in the binary dir, tell cmake to always look for includes there
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# FIND LIBRARIES ===========================================================================


# FINALIZE ==================================================================================

# generate list of source files
file(

        GLOB_RECURSE

        source_files

        src/*
)


# Set the library version in the main CMakeLists.txt
SET(GPD_MAJOR_VERSION 1)
SET(GPD_MINOR_VERSION 0)
SET(GPD_PATCH_VERSION 0)
SET(GPD_VERSION "${GPD_MAJOR_VERSION}.${GPD_MINOR_VERSION}.${GPD_PATCH_VERSION}")
SET(GPD_LIBRARY_PROPERTIES ${GPD_LIBRARY_PROPERTIES}
    VERSION "${GPD_VERSION}"
    SOVERSION "${GPD_MAJOR_VERSION}"
    SUFFIX ".so"
)


############### build the library #####################
add_library(GPD SHARED ${source_files})
target_include_directories(GPD PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_link_libraries(GPD 
         ${QT_LIBRARIES}

 	${SFML_LIBRARIES}

         ${ElementaryUtils_LIBRARIES}
 
         ${NumA++_LIBRARIES}
        
         ${PARTONS_LIBRARIES} 
MePhys Integrator Glauber ${CMAKE_SOURCE_DIR}/lib/libcuba.a)
set_target_properties(GPD PROPERTIES ${GPD_LIBRARY_PROPERTIES})
############### install the library ###################
install(TARGETS GPD DESTINATION "${CMAKE_SOURCE_DIR}/lib")

############### build the library #####################
if(STATIC_BUILD)
  add_library(GPDST STATIC ${source_files})
  target_include_directories(GPDST PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})  
  target_link_libraries(GPDST          ${QT_LIBRARIES}

  ${SFML_LIBRARIES}

  ${ElementaryUtils_LIBRARIES}

  ${NumA++_LIBRARIES}
 
  ${PARTONS_LIBRARIES} 
MePhys Integrator ${CMAKE_SOURCE_DIR}/lib/libcuba.a -Wl,-no-whole-archive)
  set_target_properties(GPDST PROPERTIES OUTPUT_NAME DIS)
  install(TARGETS GPDST DESTINATION "${CMAKE_SOURCE_DIR}/lib")
endif(STATIC_BUILD)




# define target executable
# set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# add_executable(

#         PARTONS_example

#         ${source_files}
# )

# define libraries to be linked
# target_link_libraries(

#         PARTONS_example
        
#         ${QT_LIBRARIES}

# 	${SFML_LIBRARIES}

#         ${ElementaryUtils_LIBRARIES}
 
#         ${NumA++_LIBRARIES}
        
#         ${PARTONS_LIBRARIES} 

#         Integrator 
        
#         MePhys 
# )