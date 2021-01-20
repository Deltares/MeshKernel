find_package(Doxygen)
if(NOT Doxygen_FOUND)
  message(STATUS "Doxygen not found, not building docs")
  return()
endif()

set(DOXYGEN_MAIN_PAGE ${CMAKE_CURRENT_SOURCE_DIR}/main_page.md)
set(DOXYGEN_IMAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/images)
set(DOXYGEN_INPUT_LIST
    "${PROJECT_SOURCE_DIR}/include/MeshKernel ${PROJECT_SOURCE_DIR}/include/MeshKernelApi ${DOXYGEN_MAIN_PAGE}"
)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

# This is only used for the CI
set(DOXYGEN_WARN_LOG_FILE "")

# Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

# Doxygen won't create this for us
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# meshkernel headers
file(GLOB MESHKERNEL_HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernel/*.hpp")

# meshkernelapi headers
file(GLOB MESHKERNELAPI_HEADER_LIST CONFIGURE_DEPENDS
     "${PROJECT_SOURCE_DIR}/include/MeshKernelApi/*.hpp")

# images
file(GLOB DOYGEN_IMAGES CONFIGURE_DEPENDS "${DOXYGEN_IMAGE_DIR}/*")

# Only regenerate Doxygen when the Doxyfile or public headers change
add_custom_command(
  OUTPUT ${DOXYGEN_INDEX_FILE}
  DEPENDS ${MESHKERNEL_HEADER_LIST} ${MESHKERNELAPI_HEADER_LIST}
          ${DOXYGEN_MAIN_PAGE} ${DOYGEN_IMAGES}
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
  WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/include"
  MAIN_DEPENDENCY ${DOXYFILE_OUT}
  ${DOXYFILE_IN}
  COMMENT "Generating docs")

add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})
