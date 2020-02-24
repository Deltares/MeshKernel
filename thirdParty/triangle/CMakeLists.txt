FILE( GLOB all_headers src/*.h )
source_group(include FILES ${all_headers})

FILE( GLOB all_srcs src/*.c)
source_group(source FILES ${all_srcs})

add_library(triangle 
   ${all_srcs}
   ${all_headers}
   )

install(TARGETS triangle
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib)