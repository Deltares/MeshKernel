# add version
add_subdirectory(Version)

# Add custom_memory_manager
# This target is built only if benchmarking with memory reporting is enabled
if(ENABLE_BENCHMARKING AND ENABLE_BENCHMARKING_MEM_REPORT)
  if(MEM_COLLECTION_METHOD STREQUAL "COUNT_BYTES")
    add_subdirectory(custom_memory_manager)
  elseif(MEM_COLLECTION_METHOD STREQUAL "QUERY_SYSTEM")
  add_subdirectory(memory_system_query)
  endif()
endif()

# Add test utilities
# This target is built only if unit testing is enabled
if(ENABLE_UNIT_TESTING)
  add_subdirectory(test_utils)
endif()
