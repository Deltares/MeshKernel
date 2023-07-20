#ifdef _WIN32

#pragma once

// clang-format off
#include <windows.h>
#include <psapi.h>
// clang-format on

#include "memory_monitor.hpp"

class WindowsMemoryMonitor final : public MemoryMonitor
{
public:
    WindowsMemoryMonitor()
        : m_current_process(GetCurrentProcess())
    {
    }

    ~WindowsMemoryMonitor() = default;

    // non-copyable
    WindowsMemoryMonitor(const WindowsMemoryMonitor&) = delete;
    WindowsMemoryMonitor& operator=(const WindowsMemoryMonitor&) = delete;
    WindowsMemoryMonitor(WindowsMemoryMonitor&&) = delete;
    WindowsMemoryMonitor& operator=(WindowsMemoryMonitor&&) = delete;

    uint64_t Usage(MemoryType memory_type) const override
    {
        // PROCESS_MEMORY_COUNTERS_EX pmc;
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(m_current_process,
                                 reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc),
                                 sizeof(pmc)))
        {

            switch (memory_type)
            {
            case MemoryType::Physical:
                return static_cast<uint64_t>(pmc.WorkingSetSize);
            case MemoryType::PhysicalPeak:
                return static_cast<uint64_t>(pmc.PeakWorkingSetSize);
            case MemoryType::Virtual:
                return static_cast<uint64_t>(pmc.PagefileUsage);
            case MemoryType::VirtualPeak:
                return static_cast<uint64_t>(pmc.PeakPagefileUsage);
            default:
                return 0;
            }
        }
        return 0;
    }

private:
    HANDLE m_current_process;
};

#endif