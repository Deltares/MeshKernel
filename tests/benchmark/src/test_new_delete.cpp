#include <iostream>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"
#include "custom_memory_management.hpp"
#include "test_config.hpp"
#include "test_types.hpp"

// Tests operators new and delete for single objects
TEST(NewDelete, SingleObject)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    auto pt_1 = new Point;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), size_of_point);
    delete pt_1;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), size_of_point);
    auto pt_2 = new Point;
    delete pt_2;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), 2 * size_of_point);
    auto pt_3 = new Point;
    delete pt_3;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), 3 * size_of_point);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

// Tests operators new and delete for multiple objects
TEST(NewDelete, MultipleObjects)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t const size = 5000;
    Point* const pt = new Point[size];
    delete[] pt;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), size * size_of_point);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}
