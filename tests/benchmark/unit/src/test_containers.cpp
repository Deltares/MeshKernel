#include <array>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "memory_management.hpp"
#include "test_config.hpp"
#include "test_types.hpp"

TEST(Containers, Vector)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    // let's use a scope here because we're using containers, we want them to
    // be destructed before we check for leaks towrads the end of the test
    {
        size_t const initial_size = 20;

        // construction
        std::vector<Point> vec(initial_size);

        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
        size_t const initial_bytes = vec.capacity() * size_of_point;
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), initial_bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), initial_bytes);

        // resize
        vec.resize(10); // < capacity => no realloction
        EXPECT_EQ(vec.capacity(), initial_size);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);                     // unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), initial_bytes);        // unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes); // unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), initial_bytes);       // unchaged

        // resize again
        vec.resize(200); // > capacity => realloction
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 2);
        size_t const resize_bytes = vec.capacity() * size_of_point;
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), resize_bytes);

        // clear
        vec.clear();                                                                          // capacity unchaged => no reallocation
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 2);                                    // unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);                        // unchanged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes); // unchanged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), resize_bytes);                       // unchanged

        // shrink to fit
        vec.shrink_to_fit();                                                                  // capacity becomes 0
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);                        // peak memory unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes); // unchanged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);                                  // memory is completely freed

        // push back 10 elements
        // size_t peak_bytes = CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(); // = initial_bytes + resize_bytes
        int net_heap_growth = 0;
        size_t total_allocated_bytes = CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes();
        size_t max_used_bytes = 0;
        for (int i = 0; i < 10; ++i)
        {
            size_t old_capacity = vec.capacity();
            vec.push_back(Point(i + 1, i + 2, i + 3));
            size_t new_capacity = vec.capacity();
            net_heap_growth += (int)(new_capacity - old_capacity) * size_of_point;
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), net_heap_growth);
            max_used_bytes = new_capacity * size_of_point;
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), max_used_bytes);
            if (new_capacity > old_capacity)
            {
                total_allocated_bytes += max_used_bytes;
            }
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), total_allocated_bytes);
        }

        // emplace back 20 elements
        for (int i = 0; i < 20; ++i)
        {
            size_t old_capacity = vec.capacity();
            vec.emplace_back(i - 1, i - 2, i - 3);
            size_t new_capacity = vec.capacity();
            net_heap_growth += (int)(new_capacity - old_capacity) * size_of_point;
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), net_heap_growth);
            max_used_bytes = new_capacity * size_of_point;
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), max_used_bytes);
            if (new_capacity > old_capacity)
            {
                total_allocated_bytes += max_used_bytes;
            }
            EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), total_allocated_bytes);
        }

        // remove the last element
        vec.pop_back();
        // pop_back reduces the size by 1, does not change the capacity => heap growth max bytes usedare  unchaged
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), net_heap_growth);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), total_allocated_bytes);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), max_used_bytes);
    }

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

TEST(Containers, Array)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    // let's use a scope here because we're using containers, we want them to
    // be destructed before we check for leaks towrads the end of the test
    {
        std::array<Point, 3> arr{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}};
    }

    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());

    if (config::print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}