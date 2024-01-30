
# target name undo
set(TARGET_NAME undo)

project (${TARGET_NAME} LANGUAGES CXX)

# source directory
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_LIST
  "${SOURCE_DIR}/src/SimpleMesh.cpp"
  "${SOURCE_DIR}/src/Transaction.cpp"
  "${SOURCE_DIR}/src/undo.cpp"
)

set(HEADER_LIST
  "${SOURCE_DIR}/include/AddEdgeTransaction.hpp"
  "${SOURCE_DIR}/include/AddNodeTransaction.hpp"
  "${SOURCE_DIR}/include/BaseMeshTransaction.hpp"
  "${SOURCE_DIR}/include/CompoundTransaction.hpp"
  "${SOURCE_DIR}/include/DeleteEdgeTransaction.hpp"
  "${SOURCE_DIR}/include/DeleteNodeTransaction.hpp"
  "${SOURCE_DIR}/include/MeshTypes.hpp"
  "${SOURCE_DIR}/include/SimpleMesh.hpp"
  "${SOURCE_DIR}/include/Transaction.hpp"
  "${SOURCE_DIR}/include/TransactionStack.hpp"
)

# executable
add_executable(${TARGET_NAME} ${SOURCE_LIST} ${HEADER_LIST})

set_target_properties(${TARGET_NAME} PROPERTIES LINKER_LANGUAGE CXX)

# Expose the interface of the shared lib
target_include_directories(${TARGET_NAME} PUBLIC include)
