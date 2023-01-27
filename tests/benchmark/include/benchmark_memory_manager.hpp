#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <stdexcept>
#include <stdint.h>

#include <benchmark/benchmark.h>

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
///        (implemented as a therad-safe singleton that does not register its members)
class BenchmarkMemoryManager final : public ::benchmark::MemoryManager
{
public:
    /// @brief Gets static instance of class
    static BenchmarkMemoryManager& Instance()
    {
        static BenchmarkMemoryManager instance;
        m_count++;
        if (m_count == 2)
        {
            DelayedInitialisation();
        }
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
    /// @param[ptr] The pointer to register
    /// @param[size] Size of allocated memory pointed to by pointer
    void Register(void const* const ptr, size_t size)
    {
        if (m_is_initialised)
        {
            // printf("Register\n");
            m_num_allocations++;
            m_total_allocated_bytes += size;
            m_max_bytes_used = std::max(m_max_bytes_used, m_total_allocated_bytes);
            //(*m_ptr_size_map)[ptr] = size;
            // m_ptr_size_map->insert({ptr, size});
            // m_ptr_size_map->try_emplace(ptr, size);
        }
    }

    /// @brief Unregisters an allocation
    /// @param[ptr] The pointer to deregister
    void Unregister(void const* const ptr)
    {
        if (m_is_initialised)
        {
            // printf("Unregister\n");
            m_num_deallocations++;
            // m_total_allocated_bytes -= (*m_ptr_size_map)[ptr];
            // m_ptr_size_map->erase(ptr);
        }
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
    // class is non-copyable
    BenchmarkMemoryManager(BenchmarkMemoryManager const&) = delete;
    BenchmarkMemoryManager& operator=(BenchmarkMemoryManager const&) = delete;

    static void DelayedInitialisation()
    {
        if (!m_is_initialised)
        {
            m_ptr_size_map = std::make_unique<PtrSizeMap>();
            m_is_initialised = true;
        }
    }

    int64_t m_num_allocations = 0;       ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocations = 0;     ///< The number of deallocations made in total between Start and Stop (not written to ::benchmark::MemoryManager::Result)
    int64_t m_total_allocated_bytes = 0; ///< The total memory allocated in bytes between Start and Stop
    int64_t m_max_bytes_used = 0;        ///< The peak memory use in bytes between Start and Stop
    using PtrSizeMap = std::map<void const* const, size_t>;
    inline static std::unique_ptr<PtrSizeMap> m_ptr_size_map = nullptr; ///< Map with a void pointer as key and the size of the objcet pointed to as value
    inline static size_t m_count = 0;
    inline static bool m_is_initialised = false;
};

#define MEMORY_MANAGER BenchmarkMemoryManager::Instance()
