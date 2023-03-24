# add version
add_subdirectory(Version)

# Add custom_memory_manager
# This target is built only if benchmarking with memory reporting is enabled
if(ENABLE_BENCHMARKING AND ENABLE_BENCHMARKING_MEM_REPORT)
  add_subdirectory(custom_memory_manager)
endif()

# Add test utilities
# This target is built only if unit testing is enabled
if(ENABLE_UNIT_TESTING)
  add_subdirectory(test_utils)
endif()
