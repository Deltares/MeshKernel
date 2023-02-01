#include <cstdlib>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>

#define DO_MANAGE_MEMORY true
#include "benchmark_memory_manager.hpp"
#include "custom_memory_management.hpp"

static bool constexpr print_statistics = false;

class Point
{
public:
    Point() = default;
    Point(int x, int y, int z) : m_x{x}, m_y{y}, m_z{z} {}

private:
    int m_x = -1; // 4 bytes
    int m_y = -2; // 4 bytes
    int m_z = -3; // 4 bytes
};

static size_t constexpr size_of_point = sizeof(Point); // 3 x 4 bytes = 12 bytes

TEST(MemoryManager, MallocReallocFree)
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
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);

    if (print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

// Exactly the same test as MallocReallocFree except that calloc is called in the beginning rather than malloc
TEST(MemoryManager, CallocReallocFree)
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
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);

    if (print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

// Tests operators new and delete for single objects
TEST(MemoryManager, NewDeleteSingleObject)
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
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);

    if (print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

// Tests operators new and delete for multiple objects
TEST(MemoryManager, NewDeleteMultipleObjects)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    size_t const size = 5000;
    Point* const pt = new Point[size];
    delete[] pt;
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), size * size_of_point);

    // is there a leak?
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);

    if (print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}

TEST(MemoryManager, Vector)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();
    // let's use a scope here because we're using containers, we want them to
    // be destructed before we check for leaks towrads the end of the test
    {
        size_t const initial_size = 20;

        // construction
        std::vector<Point> vec(initial_size);

        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
        size_t const initial_bytes = vec.size() * size_of_point;
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
        size_t tot_push_back_bytes = 0;
        for (int i = 0; i < 10; ++i)
        {
            vec.push_back(Point(i, i + 1, i + 2));
            // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), 1);
            // tot_push_back_bytes = vec.size() * size_of_point;
            // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), tot_push_back_bytes) << "push_back iter " << i;
            //  EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes + tot_push_back_bytes);
            // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), tot_push_back_bytes) << "push_back iter " << i;
        }

        // emplace back 20 elements
        for (int i = 0; i < 20; ++i)
        {
            vec.emplace_back(i - 1, i - 2, i - 3);
        }

        std::cout << "size " << vec.size() << " capacity  " << vec.capacity() << '\n';
        size_t const final_capacity = vec.capacity() * size_of_point;
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), final_capacity);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), resize_bytes + final_capacity);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.TotalAllocatedBytes(), initial_bytes + resize_bytes + final_capacity);
    }

    // is there a leak?
    // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.Allocations(), CUSTOM_MEMORY_MANAGER.Deallocations());
    // EXPECT_EQ(CUSTOM_MEMORY_MANAGER.NetHeapGrowth(), 0);

    if (print_statistics)
    {
        std::cout << CUSTOM_MEMORY_MANAGER;
    }
}
