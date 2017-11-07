project(progs)

# set(INCLUDE_DIRECTORIES
# ${INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
# )
message(${INCLUDE_DIRECTORIES})
include_directories( ${INCLUDE_DIRECTORIES})

set(PROG_LIST
Obstest/observables.cpp
)


foreach(it ${PROG_LIST})
    STRING(REGEX REPLACE "^.*/(.*).cpp" "\\1" bla ${it})
    message(${bla})
    set(MYPROGS ${MYPROGS} ${bla})
    add_executable(${bla} ${it} )
    target_link_libraries(${bla} ${LIBRARIES})
endforeach()

if(NOT BUILD_PROGS)
    set_target_properties(${MYPROGS} PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

if(BUILD_PROGS)
    install(TARGETS ${MYPROGS} DESTINATION "${CMAKE_SOURCE_DIR}/bin")
endif()