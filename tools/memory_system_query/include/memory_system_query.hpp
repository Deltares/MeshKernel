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
class MemorySystemQuery final : public ::benchmark::MemoryManager
{
public:
    /// @brief Gets static instance of class
    static MemorySystemQuery& Instance();

    /// @brief Starts recording allocation information
    void Start() override;

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override;

    /// @brief Gets the total memory allocated in bytes between Start and Stop
    /// @return Total memory allocated in bytes
    int64_t TotalAllocatedBytes() const;

    /// @brief Gets the peak memory use in bytes between Start and Stop
    /// @return Peak memory use in bytes between Start and Stop
    int64_t MaxBytesUsed() const;

    /// @brief Overlaod of oepartor << for printing the statistics of the class MemorySystemQuery
    /// @param[ostream] Output stream
    /// @param[ostream] Custom memory manager instance
    /// @return Output stream
    friend std::ostream& operator<<(std::ostream& ostream, MemorySystemQuery const& custom_memory_manager);

    /// @brief Gets the current statistics
    /// @return The statistics as a ormatted string
    std::string Statistics(std::string const& caller = std::string()) const;

private:
    mutable std::shared_mutex mutex;

    int64_t m_total_allocated_bytes = 0;       ///< The total memory allocated in bytes between Start and Stop
    int64_t m_max_bytes_used = TombstoneValue; ///< The peak memory use in bytes between Start and Stop

    MemorySystemQuery() = default;
    ~MemorySystemQuery() = default;
    MemorySystemQuery(MemorySystemQuery const&) = delete;
    MemorySystemQuery& operator=(MemorySystemQuery const&) = delete;

    ///
    /// The sum of the physical RAM used by the current process and the virtual RAM used by the current process.
    ///
    static int64_t GetRAMSystemUsedByCurrentProcess();

    ///
    /// On Windows, this is defined by the Working Set Size.  The working set size is defined as "The working set of a process is the set of pages in
    /// the virtual address space of the process that are currently resident in physical memory. The working set contains only pageable memory
    /// allocations; nonpageable memory allocations such as Address Windowing Extensions (AWE) or large page allocations are not included in the
    /// working set."
    /// https://msdn.microsoft.com/en-us/library/windows/desktop/ms684877(v=vs.85).aspx
    /// https://msdn.microsoft.com/en-us/library/windows/desktop/cc441804(v=vs.85).aspx
    ///
    static int64_t GetRAMPhysicalUsedByCurrentProcess();

    ///
    ///
    ///
    static int64_t GetRAMPhysicalUsedByCurrentProcessPeak();

    ///
    /// On Windows, this is defined as the commit charge.  "The Commit Charge value in bytes for this process. Commit Charge is the total amount of
    /// memory that the memory manager has committed for a running process."
    /// https://msdn.microsoft.com/en-us/library/windows/desktop/ms684877(v=vs.85).aspx
    ///
    static int64_t GetRAMVirtualUsedByCurrentProcess();
};

/// @brief Overlaod of oepartor << for printing the statistics of the class MemorySystemQuery
/// @param[ostream] Output stream
/// @param[ostream] Custom memory manager instance
/// @return Output stream
std::ostream& operator<<(std::ostream& ostream, MemorySystemQuery const& custom_memory_manager);

#define MEMORY_SYSTEM_QUERY MemorySystemQuery::Instance()
