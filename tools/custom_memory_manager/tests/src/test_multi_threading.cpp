#include "custom_memory_manager.hpp"

#include <algorithm>
#include <shared_mutex>
#include <thread>
#include <vector>

#include <gtest/gtest.h>

#include "test_config.hpp"
#include "test_types.hpp"

template <class T>
class ConcurrentVector
{
public:
    void PushBack(T const& element)
    {
        std::unique_lock lock(mutex);
        m_vector.push_back(element);
        m_bytes_used = m_vector.capacity() * sizeof(T);
    }

    size_t BytesUsed() const
    {
        std::shared_lock lock(mutex);
        return m_bytes_used;
    }

    std::vector<T> const& Get() const { return m_vector; }

private:
    mutable std::shared_mutex mutex;
    std::vector<T> m_vector;
    size_t m_bytes_used = 0;
};

TEST(Multithreading, ConcurrentVector)
{
    CUSTOM_MEMORY_MANAGER.ResetStatistics();

    {
        ConcurrentVector<int> concurrent_vector;

        CUSTOM_MEMORY_MANAGER.ResetStatistics();

        auto push_back = [&concurrent_vector](int value)
        {
            concurrent_vector.PushBack(value);
        };

        // main thread
        push_back(3);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), concurrent_vector.BytesUsed());

        // create new threads
        std::thread thread_1(push_back, 2);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), concurrent_vector.BytesUsed());
        std::thread thread_2(push_back, 0);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), concurrent_vector.BytesUsed());
        std::thread thread_3(push_back, 1);
        EXPECT_EQ(CUSTOM_MEMORY_MANAGER.MaxBytesUsed(), concurrent_vector.BytesUsed());

        thread_1.join();
        thread_2.join();
        thread_3.join();

        // sort and check if all elemenets have been inserted by the different threads
        std::vector<int> vector = concurrent_vector.Get();
        std::sort(vector.begin(), vector.end());
        EXPECT_EQ(vector, std::vector<int>({0, 1, 2, 3}));
    }
    // is there a leak?
    EXPECT_FALSE(CUSTOM_MEMORY_MANAGER.HasLeaks());
}
