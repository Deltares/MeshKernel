#include "memory_system_query.hpp"

#include <gtest/gtest.h>

#include <chrono>

#include "test_config.hpp"

TEST(Memory, AllocPhysical)
{
    const auto startUsedByCurrentProcess = MEMORY_SYSTEM_QUERY.TotalAllocatedBytes();

    // This will always be true, but the compiler won't know that, preventing the
    // allocation from happening before we want it to.
    if (std::chrono::system_clock::now() == std::chrono::time_point<std::chrono::system_clock>())
    {
        const int64_t allocAmmount = 1052672;
        volatile uint8_t* megabyte = new uint8_t[allocAmmount];

        EXPECT_NE(startUsedByCurrentProcess, MEMORY_SYSTEM_QUERY.TotalAllocatedBytes());

        EXPECT_LT(startUsedByCurrentProcess + allocAmmount, MEMORY_SYSTEM_QUERY.TotalAllocatedBytes())
            << "We did an allocation, but are using less RAM. Start: " << startUsedByCurrentProcess << ", Alloc: " << allocAmmount;

        // Assume we are not swapping out to disk
        ASSERT_GT(MEMORY_SYSTEM_QUERY.TotalAllocatedBytes(), startUsedByCurrentProcess);

        const auto memoryDelta = MEMORY_SYSTEM_QUERY.TotalAllocatedBytes() - startUsedByCurrentProcess;
        EXPECT_NEAR(static_cast<double>(MEMORY_SYSTEM_QUERY.TotalAllocatedBytes() - startUsedByCurrentProcess),
                    static_cast<double>(memoryDelta),
                    4096.0)
            << "Memory Delta: " << memoryDelta;

        delete[] megabyte;
    }
}
