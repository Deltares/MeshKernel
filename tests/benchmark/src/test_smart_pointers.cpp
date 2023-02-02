#include <iostream>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"
#include "custom_memory_management.hpp"
#include "test_config.hpp"
#include "test_types.hpp"

// Tests operators new and delete for single objects
TEST(SmartPointer, SharedSingleObject)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    // let's use a scope here because we're using containers, we want them to
    // be destructed before we check for leaks towrads the end of the test
    {
        std::shared_ptr<Point> ptr = std::make_shared<Point>(-6, -6, -6);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
        long constexpr size_of_ref_count = sizeof(long);
        size_t const bytes = size_of_point + sizeof(ptr) + size_of_ref_count;
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), bytes);
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Deallocations(), 1);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

TEST(SmartPointer, UniqueMultipleObject)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    // let's use a scope here because we're using containers, we want them to
    // be destructed before we check for leaks towrads the end of the test
    {
        size_t const size = 3;
        std::unique_ptr<Point[]> ptr = std::make_unique<Point[]>(size);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
        size_t const bytes = 3 * size_of_point;
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), bytes);
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Deallocations(), 1);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}
