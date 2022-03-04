# Set library files
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")
 
file(GLOB SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/MeshKernel/*.cpp")

file(GLOB AVERAGING_STRATEGIES_HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/MeshKernel/AveragingStrategies/*.hpp")
 
file(GLOB AVERAGING_STRATEGIES_SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/MeshKernel/AveragingStrategies/*.cpp")

file(GLOB CURVILINEAR_GRID_HEADER_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/include/MeshKernel/CurvilinearGrid/*.hpp")

file(GLOB CURVILINEAR_GRID_SOURCE_LIST CONFIGURE_DEPENDS "${PROJECT_SOURCE_DIR}/src/MeshKernel/CurvilinearGrid/*.cpp")

# Create the static lib
add_library(MeshKernel STATIC ${SOURCE_LIST} ${HEADER_LIST} 
                              ${AVERAGING_STRATEGIES_HEADER_LIST} ${AVERAGING_STRATEGIES_SOURCE_LIST}
                              ${CURVILINEAR_GRID_HEADER_LIST} ${CURVILINEAR_GRID_SOURCE_LIST})

# Make sure that coverage information is produced when using gcc
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set_target_properties(
    MeshKernel
    PROPERTIES COMPILE_FLAGS
               "-fprofile-arcs -ftest-coverage -static-libgcc -static-libstdc++"
  )
  set_target_properties(MeshKernel PROPERTIES LINK_FLAGS
                                              "-ldl -lgcov --coverage")
endif()

# Expose the interface of the static lib
target_include_directories(MeshKernel PUBLIC "${PROJECT_SOURCE_DIR}/include")

# Add target link dependencies
target_link_libraries(MeshKernel LINK_PUBLIC ${Boost_LIBRARIES} triangle
                      OpenMP::OpenMP_CXX)

# Add precompiled header to speed up compilation
target_precompile_headers(MeshKernel PRIVATE <map>
                                             <string> 
                                             <unordered_map> 
                                             <vector> 
                                             <memory> 
                                             <algorithm>
                                             <cmath>
                                             <numeric>
                                             <tuple>
                                             <limits>
                                             <utility>
                                             <array>
                                             <initializer_list>
                                             <stdexcept> 
                                             <cassert>)

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernel"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
  
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernel/AveragingStrategies"
  PREFIX "Header Files/AveragingStrategies"
  FILES ${AVERAGING_STRATEGIES_HEADER_LIST})
  
source_group(
  TREE "${PROJECT_SOURCE_DIR}/src/MeshKernel/AveragingStrategies"
  PREFIX "Source Files/AveragingStrategies"
  FILES ${AVERAGING_STRATEGIES_SOURCE_LIST})
  
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include/MeshKernel/CurvilinearGrid"
  PREFIX "Header Files/CurvilinearGrid"
  FILES ${CURVILINEAR_GRID_HEADER_LIST})
  
source_group(
  TREE "${PROJECT_SOURCE_DIR}/src/MeshKernel/CurvilinearGrid"
  PREFIX "Source Files/CurvilinearGrid"
  FILES ${CURVILINEAR_GRID_SOURCE_LIST})

install(TARGETS MeshKernel)