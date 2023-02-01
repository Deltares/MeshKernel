#pragma once

#include <algorithm>
#include <iostream>
#include <stdint.h>

#include <benchmark/benchmark.h>

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
///        (implemented as a therad-safe singleton)
class CustomMemoryManager final : public ::benchmark::MemoryManager
{
public:
    /// @brief Gets static instance of class
    static CustomMemoryManager& Instance()
    {
        static CustomMemoryManager instance;
        return instance;
    }

    /// @brief Starts recording allocation information
    void Start() override
    {
        ResetStatistics();
    }

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override
    {
        result->num_allocs = m_num_allocations;
        result->max_bytes_used = m_max_bytes_used;
        result->total_allocated_bytes = m_total_allocated_bytes;
        result->net_heap_growth = m_net_heap_growth;
    }

    /// @brief Registers an allocation
    /// @param[size] Size of allocated memory pointed to by pointer
    void Register(int64_t size, bool incrrement_num_allocations = true)
    {
        if (incrrement_num_allocations)
        {
            m_num_allocations++;
        }
        m_total_allocated_bytes += size;
        m_net_heap_growth += size;
        m_max_bytes_used = std::min(m_max_bytes_used, m_net_heap_growth);
    }

    /// @brief Unregisters an allocation
    /// @param[size] Size of allocated memory pointed to by pointer
    void Unregister(int64_t size)
    {
        m_num_deallocations++;
        m_net_heap_growth -= size;
        m_max_bytes_used = std::max(m_max_bytes_used, m_net_heap_growth);
    }

    /// @brief Gets the total number of allocations
    /// @return Number of allocations
    int64_t Allocations() const { return m_num_allocations; }

    /// @brief Gets the total number of deallocations
    /// @return Number of deallocations
    int64_t Deallocations() const { return m_num_deallocations; }

    /// @brief Gets the total memory allocated in bytes between Start and Stop
    /// @return Total memory allocated in bytes
    int64_t TotalAllocatedBytes() const { return m_total_allocated_bytes; }

    /// @brief Gets the peak memory use in bytes between Start and Stop
    /// @return Peak memory use in bytes between Start and Stop
    int64_t MaxBytesUsed() const { return m_max_bytes_used; }

    /// @brief Gets the net changes in memory in bytes between Start and Stop
    /// @return Net changes in memory in bytes between Start and Stop
    int64_t NetHeapGrowth() const { return m_net_heap_growth; }

    /// @brief Resets the memory statistics
    void ResetStatistics()
    {
        m_num_allocations = 0;
        m_num_deallocations = 0;
        m_total_allocated_bytes = 0;
        m_max_bytes_used = TombstoneValue;
        m_net_heap_growth = 0;
    }

    friend std::ostream& operator<<(std::ostream& ostream, CustomMemoryManager const& custom_memory_manager);

private:
    CustomMemoryManager() = default;
    ~CustomMemoryManager() = default;
    CustomMemoryManager(CustomMemoryManager const&) = delete;
    CustomMemoryManager& operator=(CustomMemoryManager const&) = delete;

    int64_t m_num_allocations = 0;             ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocations = 0;           ///< The number of deallocations made in total between Start and Stop (not written to ::benchmark::MemoryManager::Result)
    int64_t m_total_allocated_bytes = 0;       ///< The total memory allocated in bytes between Start and Stop
    int64_t m_max_bytes_used = TombstoneValue; ///< The peak memory use in bytes between Start and Stop
    int64_t m_net_heap_growth = 0;             ///< The net changes in memory in bytes between Start and Stop
};

inline static std::ostream& operator<<(std::ostream& ostream, CustomMemoryManager const& custom_memory_manager)
{
    ostream << "Current memory manager statistics:"
            << "\nNumber of allocations  : " << custom_memory_manager.m_num_allocations
            << "\nNumber of deallocations: " << custom_memory_manager.m_num_deallocations
            << "\nTotal allocated bytes  : " << custom_memory_manager.m_total_allocated_bytes
            << "\nMax bytes used         : " << custom_memory_manager.m_max_bytes_used
            << "\nNet heap growth        : " << custom_memory_manager.m_net_heap_growth
            << '\n';
    return ostream;
}

#define CUSTOM_MEMORY_MANAGER CustomMemoryManager::Instance()
