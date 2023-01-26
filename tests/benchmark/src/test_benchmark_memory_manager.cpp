#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>

#include <crtdbg.h>

#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"

#define DO_MANAGE_MEMORY true
#include "custom_memory_management.cpp"

class MemoryLeakDetector
{
public:
    MemoryLeakDetector()
    {
        _CrtMemCheckpoint(&memState_);
    }

    ~MemoryLeakDetector()
    {
        _CrtMemState stateNow, stateDiff;
        _CrtMemCheckpoint(&stateNow);
        int diffResult = _CrtMemDifference(&stateDiff, &memState_, &stateNow);
        if (diffResult)
        {
            reportFailure((int)stateDiff.lSizes[1]);
        }
    }

private:
    void reportFailure(unsigned int unfreedBytes)
    {
        FAIL() << " >>>>>>>>> Memory leak of " << unfreedBytes << " byte(s) detected.";
    }
    _CrtMemState memState_;
};

struct Point
{
    int x = -1;
    int y = -1;
};

TEST(MemoryManager, NewDeleteSingleObject)
{
    MEMORY_MANAGER.ResetAllocations();   // remove me
    MEMORY_MANAGER.ResetDeallocations(); // remove me
    MemoryLeakDetector leakDetector;
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
    MemoryLeakDetector leakDetector;
    Point* const pt = new Point[5];
    delete[] pt;
    // EXPECT_EQ(MEMORY_MANAGER.Allocations(), MEMORY_MANAGER.Deallocations());
    std::cout << "TRACKMEM: " << MEMORY_MANAGER.Allocations() << ' ' << MEMORY_MANAGER.Deallocations() << '\n';
}
*/

#undef DO_MANAGE_MEMORY
