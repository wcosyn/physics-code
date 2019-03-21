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

# find libraries: Qt4
find_package(Qt4 REQUIRED QtCore QtSql)

# add some useful macros and variables
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

# find libraries: additional modules
message(${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

# find libraries: SFML
find_package(SFML COMPONENTS system REQUIRED)

# find libraries: ElementaryUtils
find_package(ElementaryUtils REQUIRED)

# find libraries: NumA++
find_package(NumA++ REQUIRED)

# find libraries: PARTONS
find_package(PARTONS REQUIRED)

# directories containing headers 
include_directories(${QT_INCLUDE_DIRS} ${SFML_INCLUDE_DIR} ${ELEMENTARY_UTILS_INCLUDE_DIR} ${NUMA_INCLUDE_DIR} ${PARTONS_INCLUDE_DIR})

# FINALIZE ==================================================================================

# generate list of source files
file(

        GLOB_RECURSE

        source_files

        src/*
)

# define target executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(

        PARTONS_example

        ${source_files}
)

# define libraries to be linked
target_link_libraries(

        PARTONS_example
        
        ${QT_LIBRARIES}

	${SFML_LIBRARIES}

        ${ELEMENTARY_UTILS_LIBRARIES}
 
        ${NUMA_LIBRARIES}
        
        ${PARTONS_LIBRARIES} 

        Integrator 
        
        MePhys 
)