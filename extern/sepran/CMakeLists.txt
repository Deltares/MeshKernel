enable_language(Fortran)

set(target_name Sepran)

add_library(${target_name} STATIC)

set(
  TARGET_SRC_LIST
  chsort.for
  mshconstants.f90
  mshdummymethods.f90
  msh401.for
  msh402.for
  msh403.for
  msh406.for
  msh416.for
  mshchkstapl.for
  mshcopyboun.for
  mshcrossline.for
  mshcrossline1.for
  mshcurvinters.for
  mshcurvinters1.for
  mshcurvinters2.for
  msho01.for
  msho02.for
  msho03.for
  msho04.for
  msho05.for
  msho06.for
  msho07.for
  msho08.for
  msho09.for
  msho10.for
  msho11.for
  msho12.for
  msho13.for
  msho14.for
  msho15.for
  msho16.for
  msho17.for
  msho18.for
  msho19.for
  msho20.for
  msho21.for
  msho22.for
  msho24.for
  msho25.for
  msho26.for
  msho27.for
  msho28.for
  msho29.for
  msho2d.for
  msho30.for
  msho31.for
  msho32.for
  msho33.for
  msho34.for
  msho35.for
  msho36.for
  msho38.for
  msho39.for
  msho40.for
  msho41.for
  msho42.for
  msho75.for
  mshoce.for
  mshtrans2dsur.for
)

target_sources(${target_name} PRIVATE ${TARGET_SRC_LIST})

set_target_properties(
  ${target_name}
  PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(${target_name} PRIVATE -w)
  target_link_libraries(${target_name} INTERFACE gfortran)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(${target_name} PRIVATE -w)
endif()

source_group("Source Files" FILES ${TARGET_SRC_LIST})