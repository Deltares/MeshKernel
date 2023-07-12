#pragma once

#include <cstdint>

enum class MemoryType
{
    Physical = 0,
    PhysicalPeak = 1,
    Virtual = 2,
    VirtualPeak = 3
};

class MemoryMonitor
{
public:
    MemoryMonitor() = default;
    virtual ~MemoryMonitor() = default;
    // non-copyable
    MemoryMonitor(const MemoryMonitor& other) = delete;
    MemoryMonitor& operator=(const MemoryMonitor& other) = delete;
    MemoryMonitor(MemoryMonitor&& other) = delete;
    MemoryMonitor& operator=(MemoryMonitor&& other) = delete;

    /// @brief Gets the memory usage my type
    /// @param memory_type The type of requested memory usage
    /// @return            The memory usage
    virtual uint64_t Usage(MemoryType memory_type) const = 0;

    /// @brief Convenience method that returns the total current physical and virtual usage
    /// @return The total current physical and virtual usage
    uint64_t CurrentUsage() const
    {
        return Usage(MemoryType::Physical) + Usage(MemoryType::Virtual);
    }

    /// @brief Convenience method that returns the total peak physical and virtual usage
    /// @return The total peak physical and virtual usage
    uint64_t PeakUsage() const
    {
        return Usage(MemoryType::PhysicalPeak) + Usage(MemoryType::VirtualPeak);
    }

private:
    inline static uint64_t constexpr kilobyte = 1024;
};
