find_package(Doxygen)
if(NOT Doxygen_FOUND)
  message(STATUS "Doxygen not found, not building docs")
  return()
endif()

find_package(Sphinx)
if(NOT Sphinx_FOUND)
  message(STATUS "Sphinx not found, not building docs")
  return()
endif()

set(DOXYGEN_INPUT_DIR "${PROJECT_SOURCE_DIR}/include/MeshKernel")
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) # Doxygen won't create this for us
add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE}
  DEPENDS "${HEADER_LIST}"
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/include"
  MAIN_DEPENDENCY ${DOXYFILE_OUT}
  ${DOXYFILE_IN}
  COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
