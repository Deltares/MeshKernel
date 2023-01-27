#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"

struct Point
{
    int x = -1;
    int y = -1;
};

TEST(MemoryManager, NewDeleteSingleObject)
{
    MEMORY_MANAGER.ResetAllocations();   // remove me
    MEMORY_MANAGER.ResetDeallocations(); // remove me
    auto pt = new Point;
    delete pt;
    // std::vector<int> vec(100, 10);
    //  vec.push_back(1);
    // vec.resize(30);
    // vec.resize(3);
    //  vec.clear();
    //  vec.shrink_to_fit();
    //   EXPECT_EQ(MEMORY_MANAGER.Allocations(), MEMORY_MANAGER.Deallocations());
    std::cout << "MemoryManager.NewDeleteSingleObject TRACKMEM: "
              << MEMORY_MANAGER.Allocations() << ' '
              << MEMORY_MANAGER.Deallocations() << '\n';
}

/*
TEST(MemoryManager, NewDeleteMultipleObjects)
{
    Point* const pt = new Point[5];
    delete[] pt;
    // EXPECT_EQ(MEMORY_MANAGER.Allocations(), MEMORY_MANAGER.Deallocations());
    std::cout << "TRACKMEM: " << MEMORY_MANAGER.Allocations() << ' ' << MEMORY_MANAGER.Deallocations() << '\n';
}
*/

#undef DO_MANAGE_MEMORY
