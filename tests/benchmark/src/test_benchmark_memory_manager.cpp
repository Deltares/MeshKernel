#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"
// #define DO_MANAGE_MEMORY true
#include "custom_memory_management.hpp"

BenchmarkMemoryManager const& memory_manager = MEMORY_MANAGER;

struct Point
{
    int x = -1;
    int y = -1;
};

TEST(MemoryManager, NewDeleteSingleObject)
{
    Point* const pt = new Point;
    delete pt;
    EXPECT_EQ(memory_manager.Allocations(), memory_manager.Deallocations());
    std::cout << "TRACKMEM: " << memory_manager.Allocations() << ' ' << memory_manager.Deallocations() << '\n';
}

TEST(MemoryManager, NewDeleteMultipleObjects)
{
    Point* const pt = new Point[5];
    delete[] pt;
    EXPECT_EQ(memory_manager.Allocations(), memory_manager.Deallocations());
    std::cout << "TRACKMEM: " << memory_manager.Allocations() << ' ' << memory_manager.Deallocations() << '\n';
}
