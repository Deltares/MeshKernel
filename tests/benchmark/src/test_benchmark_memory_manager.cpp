#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#define DO_MANAGE_MEMORY true
#include "benchmark_memory_manager.hpp"
#include "custom_memory_management.hpp"

struct Point
{
    int x = -1;
    int y = -2;
    int z = -3;
};

TEST(MemoryManager, XAlloc)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    Point* pt = static_cast<Point*>(calloc(1, sizeof(struct Point)));
    *pt = Point();
    ASSERT_TRUE(pt != nullptr);
    EXPECT_EQ(pt->x, -1);
    EXPECT_EQ(pt->y, -2);
    EXPECT_EQ(pt->z, -3);
    std::cout << pt << '\n';
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
    void* new_pt = realloc(pt, 10);
    ASSERT_TRUE(new_pt != nullptr);
    pt = static_cast<Point*>(new_pt);
    ASSERT_TRUE(pt != nullptr);
    std::cout << pt << '\n';
    // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 2);
    free(pt);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Deallocations(), 1);
}

TEST(MemoryManager, NewDeleteSingleObject)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    auto pt_1 = new Point;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), sizeof(Point));
    delete pt_1;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), sizeof(Point));
    auto pt_2 = new Point;
    delete pt_2;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), 2 * sizeof(Point));
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
}

TEST(MemoryManager, NewDeleteMultipleObjects)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t const size = 5000;
    Point* const pt = new Point[size];
    delete[] pt;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), size * sizeof(Point));
}

TEST(MemoryManager, Vector)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t const initial_size = 20;
    std::vector<Point> vec(initial_size);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
    size_t const initial_bytes = vec.size() * sizeof(Point);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), initial_bytes);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), initial_bytes);
    //
    vec.resize(10); // < capacity => no realloction
    EXPECT_EQ(vec.capacity(), initial_size);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);                     // unchaged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), initial_bytes);        // unchaged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes); // unchaged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), initial_bytes);       // unchaged
    //
    vec.resize(200); // > capacity => realloction
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 2);
    size_t const resize_bytes = vec.size() * sizeof(Point);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), resize_bytes);
    //
    vec.clear();                                                                          // capacity unchaged => no reallocation
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes); // unchanged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), resize_bytes);                       // unchanged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);                        // unchanged
    //
    vec.shrink_to_fit(); // capacity becomes 0
    // size_t const clear_bytes = vec.size() * sizeof(Point);
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes); // unchanged
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);                                  // memory is completely freed
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes);                        // peak memory unchaged
}
