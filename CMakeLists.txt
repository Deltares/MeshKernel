project( nodeAdminstration )

cmake_minimum_required(VERSION 3.1)

FIND_PACKAGE( Boost )
FIND_PACKAGE( Cgal )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

ADD_EXECUTABLE( nodeAdministration nodeAdministration.cpp )
TARGET_LINK_LIBRARIES( nodeAdministration LINK_PUBLIC ${Boost_LIBRARIES} CGAL::CGAL )