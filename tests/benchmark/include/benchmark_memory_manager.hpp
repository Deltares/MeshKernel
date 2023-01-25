#pragma once

#include <benchmark/benchmark.h>

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
class BenchmarkMemoryManager final : public benchmark::MemoryManager
{
public:
    /// @brief Gets static instance of class
    static BenchmarkMemoryManager& Instance()
    {
        static BenchmarkMemoryManager instance;
        return instance;
    }

    /// @brief Starts recording allocation information
    void Start() override
    {
        m_num_allocs = 0;
        m_num_deallocs = 0;
        m_total_allocated_bytes = 0;
        m_max_bytes_used = 0;
    }

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override
    {
        result->num_allocs = m_num_allocs;
        result->total_allocated_bytes = m_total_allocated_bytes;
        result->max_bytes_used = m_max_bytes_used;
    }

    /// @brief Registers an allocation
    /// @param[size] Size of allocated memory
    void RegisterAllocation(size_t size)
    {
        m_num_allocs++;
        m_total_allocated_bytes += size;
        m_max_bytes_used = std::max(m_max_bytes_used, m_total_allocated_bytes);
    }

    /// @brief Registers an deallocation
    /// @param[ptr] size of allocated memory
    void RegisterDeallocation(void const* const ptr)
    {
        m_num_deallocs++;
        m_total_allocated_bytes -= sizeof(ptr); // this is incorrect!!!
    }

    /// @brief Gets total number of allocations
    /// @return Number of allocations
    size_t Allocations() const { return m_num_allocs; }

    /// @brief Gets total number of deallocations
    /// @return Number of deallocations
    size_t Deallocations() const { return m_num_deallocs; }

private:
    BenchmarkMemoryManager() = default;
    ~BenchmarkMemoryManager() = default;
    BenchmarkMemoryManager(BenchmarkMemoryManager const&) = delete;
    BenchmarkMemoryManager& operator=(BenchmarkMemoryManager const&) = delete;

    int64_t m_num_allocs = 0;            ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocs = 0;          ///< The number of deallocations made in total between Start and Stop
    int64_t m_total_allocated_bytes = 0; ///< The total memory allocated, in bytes, between Start and Stop
    int64_t m_max_bytes_used = 0;        ///< The peak memory use between Start and Stop
};

#define MEMORY_MANAGER BenchmarkMemoryManager::Instance()
