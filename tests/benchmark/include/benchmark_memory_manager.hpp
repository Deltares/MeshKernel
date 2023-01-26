#pragma once

#include <map>
#include <stdint.h>

#include <benchmark/benchmark.h>

// #include "custom_memory_management.hpp"

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
///        (implemented as a therad-safe singleton)
class BenchmarkMemoryManager final : public ::benchmark::MemoryManager
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
        m_num_allocations = 0;
        m_num_deallocations = 0;
        m_total_allocated_bytes = 0;
        m_max_bytes_used = 0;
    }

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override
    {
        result->num_allocs = m_num_allocations;
        result->total_allocated_bytes = m_total_allocated_bytes;
        result->max_bytes_used = m_max_bytes_used;
    }

    /// @brief Registers an allocation
    /// @param[size] Size of allocated memory
    void RegisterAllocation(void const* const ptr, size_t size)
    {
        m_num_allocations++;
        m_total_allocated_bytes += size;
        m_max_bytes_used = std::max(m_max_bytes_used, m_total_allocated_bytes);
        // m_sizeof.insert({ptr, size});
    }

    /// @brief Registers an deallocation
    /// @param[ptr] size of allocated memory
    void RegisterDeallocation(void const* const ptr)
    {
        m_num_deallocations++;
        m_total_allocated_bytes -= 0; // m_sizeof[ptr];
        // m_sizeof.erase(ptr);
    }

    /// @brief Gets the total number of allocations
    /// @return Number of allocations
    size_t Allocations() const { return m_num_allocations; }

    /// @brief Gets the total number of deallocations
    /// @return Number of deallocations
    size_t Deallocations() const { return m_num_deallocations; }

    /// @brief Gets the total memory allocated in bytes between Start and Stop
    /// @return Total memory allocated in bytes
    size_t TotalAllocatedBytes() const { return m_total_allocated_bytes; }

    /// @brief Gets the peak memory use in bytes between Start and Stop
    /// @return Peak memory use in bytes between Start and Stop
    size_t MaxBytesUsed() const { return m_max_bytes_used; }

    // remove me: used to overcome gtest leaks
    void ResetAllocations() { m_num_allocations = 0; }

    // remove me: used to overcome gtest leaks
    void ResetDeallocations() { m_num_deallocations = 0; }

private:
    BenchmarkMemoryManager() = default;
    ~BenchmarkMemoryManager() = default;
    BenchmarkMemoryManager(BenchmarkMemoryManager const&) = delete;
    BenchmarkMemoryManager& operator=(BenchmarkMemoryManager const&) = delete;

    int64_t m_num_allocations = 0;       ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocations = 0;     ///< The number of deallocations made in total between Start and Stop (not written to ::benchmark::MemoryManager::Result)
    int64_t m_total_allocated_bytes = 0; ///< The total memory allocated in bytes between Start and Stop
    int64_t m_max_bytes_used = 0;        ///< The peak memory use in bytes between Start and Stop
    // std::map<void const* const, size_t> m_sizeof; ///< Map with a void pointer as key and the size of the objcet pointed to as value
};

#define MEMORY_MANAGER BenchmarkMemoryManager::Instance()
