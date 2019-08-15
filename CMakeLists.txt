project( GridGeom )
cmake_minimum_required(VERSION 3.1)

FIND_PACKAGE( Boost )
FIND_PACKAGE( Cgal )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

FILE( GLOB all_headers src/*.h )
source_group(include FILES ${all_headers})

FILE( GLOB all_srcs src/*.cpp)
source_group(source FILES ${all_srcs})

add_executable( GridGeom
  ${all_srcs}
  ${all_headers}
)

target_link_libraries(GridGeom LINK_PUBLIC ${Boost_LIBRARIES} CGAL::CGAL)