set(TRIANGLE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/triangle.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/TriangleFacade.c
    ${CMAKE_CURRENT_SOURCE_DIR}/src/triangle.h)

add_definitions(-DTRILIBRARY)
add_definitions(-DNO_TIMER)
add_definitions(-w)

add_library(triangle STATIC ${TRIANGLE_SOURCES})

set(TRIANGLE_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(TRIANGLE_LIBRARIES triangle)
