#include <iostream>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"
#include "custom_memory_management.hpp"
#include "test_config.hpp"
#include "test_types.hpp"

TEST(CAllocators, MallocReallocFree)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t num_objects = 10;
    size_t tot_num_objects = num_objects;
    Point* ptr = static_cast<Point*>(malloc(num_objects * size_of_point));
    ASSERT_TRUE(ptr != nullptr);
    size_t allocs = CUSTOM_MEMORY_MANAGER.Allocations();
    EXPECT_EQ(allocs, 1);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // try to expand, be conservative how much more mem is requested
    num_objects = 11;
    void* new_ptr_1 = realloc(ptr, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_1 != nullptr);
    if (new_ptr_1 != ptr)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 10;
        EXPECT_EQ(new_ptr_1, ptr); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // try to expand again, be conservative how much more mem is requested
    num_objects = 12;
    void* new_ptr_2 = realloc(new_ptr_1, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_2 != nullptr);
    if (new_ptr_2 != new_ptr_1)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 11;
        EXPECT_EQ(new_ptr_2, new_ptr_1); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // shrink
    num_objects = 8;
    tot_num_objects += num_objects - 12;
    void* new_ptr_3 = realloc(new_ptr_2, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_3 != nullptr);
    EXPECT_EQ(new_ptr_3, new_ptr_2);                        // shrinking does not change the address
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs); // note allocs was not incremented

    // try to reallocate, request a lot more memory
    num_objects = 100000;
    void* new_ptr_4 = realloc(new_ptr_3, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_4 != nullptr);
    if (new_ptr_4 != new_ptr_3)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 12;
        EXPECT_EQ(new_ptr_4, new_ptr_3); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    free(new_ptr_4);

    // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

// Exactly the same test as MallocReallocFree except that calloc is called in the beginning rather than malloc
TEST(CAllocators, CallocReallocFree)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t num_objects = 10;
    size_t tot_num_objects = num_objects;
    Point* ptr = static_cast<Point*>(calloc(num_objects, size_of_point));
    ASSERT_TRUE(ptr != nullptr);
    size_t allocs = CUSTOM_MEMORY_MANAGER.Allocations();
    EXPECT_EQ(allocs, 1);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // try to expand, be conservative how much more mem is requested
    num_objects = 11;
    void* new_ptr_1 = realloc(ptr, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_1 != nullptr);
    if (new_ptr_1 != ptr)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 10;
        EXPECT_EQ(new_ptr_1, ptr); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // try to expand again, be conservative how much more mem is requested
    num_objects = 12;
    void* new_ptr_2 = realloc(new_ptr_1, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_2 != nullptr);
    if (new_ptr_2 != new_ptr_1)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 11;
        EXPECT_EQ(new_ptr_2, new_ptr_1); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // shrink
    num_objects = 8;
    tot_num_objects += num_objects - 12;
    void* new_ptr_3 = realloc(new_ptr_2, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_3 != nullptr);
    EXPECT_EQ(new_ptr_3, new_ptr_2);                        // shrinking does not change the address
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs); // note allocs was not incremented

    // try to reallocate, request a lot more memory
    num_objects = 100000;
    void* new_ptr_4 = realloc(new_ptr_3, num_objects * size_of_point);
    ASSERT_TRUE(new_ptr_4 != nullptr);
    if (new_ptr_4 != new_ptr_3)
    {
        allocs++;
        tot_num_objects += num_objects;
    }
    else
    {
        tot_num_objects += num_objects - 12;
        EXPECT_EQ(new_ptr_4, new_ptr_3); // expansion does not change the address
    }
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), allocs);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    free(new_ptr_4);

    // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), tot_num_objects * size_of_point);

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}
