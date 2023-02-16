#pragma once

#include <atomic>
#include <mutex>
#include <ostream>
#include <shared_mutex>
#include <string>

#include <stdint.h>

#include <benchmark/benchmark.h>

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
///        (implemented as a therad-safe singleton)
class CustomMemoryManager final : public ::benchmark::MemoryManager
{
public:
    /// @brief Gets static instance of class
    static CustomMemoryManager& Instance();

    /// @brief Starts recording allocation information
    void Start() override;

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override;

    /// @brief Custom std::malloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
    /// @param[size] Number of bytes of uninitialized storage to be allocated
    /// @return Pointer to the beginning of newly allocated memory
    void* Alloc(size_t size);

    /// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
    ///        which registers the allocation in a global BenchmarkMemoryManager object
    /// @param[size] Number of bytes of uninitialized storage to be allocated
    /// @param[alignment] Specifies the alignment
    /// @return The pointer to the beginning of newly allocated memory
    void* AlignedAlloc(size_t size, size_t alignment);

    /// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
    /// @param[ptr] Pointer to the memory block to deallocate
    void Free(void* ptr);

    /// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
    ///        and registers the allocation in a global BenchmarkMemoryManager object
    /// @param[ptr] Pointer to the memory block to deallocate
    void AlignedFree(void* ptr);

    /// @brief Gets the total number of allocations
    /// @return Number of allocations
    int64_t Allocations() const;

    /// @brief Gets the total number of deallocations
    /// @return Number of deallocations
    int64_t Deallocations() const;

    /// @brief Gets the total memory allocated in bytes between Start and Stop
    /// @return Total memory allocated in bytes
    int64_t TotalAllocatedBytes() const;

    /// @brief Gets the peak memory use in bytes between Start and Stop
    /// @return Peak memory use in bytes between Start and Stop
    int64_t MaxBytesUsed() const;

    /// @brief Gets the net changes in memory in bytes between Start and Stop
    /// @return Net changes in memory in bytes between Start and Stop
    int64_t NetHeapGrowth() const;

    /// @brief Checks if there are leaks based on the net between the number of allocations
    ///        and deallocations or the net heap growth
    /// @return Boolean indicating whether a leak has been detected
    bool HasLeaks() const;

    /// @brief Resets the memory statistics
    void ResetStatistics();

    /// @brief Overlaod of oepartor << for printing the statistics of the class CustomMemoryManager
    /// @param[ostream] Output stream
    /// @param[ostream] Custom memory manager instance
    /// @return Output stream
    friend std::ostream& operator<<(std::ostream& ostream, CustomMemoryManager const& custom_memory_manager);

    /// @brief Gets the current statistics
    /// @return The statistics as a ormatted string
    std::string Statistics(std::string const& caller = std::string()) const;

private:
    mutable std::shared_mutex mutex;

    int64_t m_num_allocations = 0;             ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocations = 0;           ///< The number of deallocations made in total between Start and Stop (not written to ::benchmark::MemoryManager::Result)
    int64_t m_total_allocated_bytes = 0;       ///< The total memory allocated in bytes between Start and Stop
    int64_t m_max_bytes_used = TombstoneValue; ///< The peak memory use in bytes between Start and Stop
    int64_t m_net_heap_growth = 0;             ///< The net changes in memory in bytes between Start and Stop

    CustomMemoryManager() = default;
    ~CustomMemoryManager() = default;
    CustomMemoryManager(CustomMemoryManager const&) = delete;
    CustomMemoryManager& operator=(CustomMemoryManager const&) = delete;

    /// @brief Registers an allocation
    /// @param[size] Size of allocated memory pointed to by pointer
    void Register(int64_t size);

    /// @brief Unregisters an allocation
    /// @param[size] Size of allocated memory pointed to by pointer
    void Unregister(int64_t size);

    /// @brief Returns the size, in bytes, of a memory block allocated in the heap
    /// @param[ptr] Pointer to the memory block
    size_t MemoryBlockSize(void* ptr);
};

/// @brief Overlaod of oepartor << for printing the statistics of the class CustomMemoryManager
/// @param[ostream] Output stream
/// @param[ostream] Custom memory manager instance
/// @return Output stream
std::ostream& operator<<(std::ostream& ostream, CustomMemoryManager const& custom_memory_manager);

#define CUSTOM_MEMORY_MANAGER CustomMemoryManager::Instance()
