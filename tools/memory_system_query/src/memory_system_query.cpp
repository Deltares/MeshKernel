#include "memory_system_query.hpp"

#include <map>
#include <sstream>
#include <string>
#include <utility>

#if defined(_WIN32)
// clang-format off
#include <windows.h>
#include <psapi.h>
// clang-format on
#elif defined(__linux__)
#include "memory_monitor_linux.hpp"
#else
#error "Unsupported platform"
#endif

MemorySystemQuery& MemorySystemQuery::Instance()
{
    static MemorySystemQuery instance;
    return instance;
}

void MemorySystemQuery::Start()
{
    std::unique_lock lock(mutex);
    m_total_allocated_bytes = CurrentMemoryUsage();
    m_max_bytes_used = PeakMemoryUsage();
}

void MemorySystemQuery::Stop(Result& result)
{
    std::unique_lock lock(mutex);
    result.total_allocated_bytes = CurrentMemoryUsage() - m_total_allocated_bytes;
    result.max_bytes_used = PeakMemoryUsage() - m_max_bytes_used;
}

int64_t MemorySystemQuery::TotalAllocatedBytes() const
{
    std::shared_lock lock(mutex);
    return CurrentMemoryUsage() - m_total_allocated_bytes;
}

int64_t MemorySystemQuery::MaxBytesUsed() const
{
    std::shared_lock lock(mutex);
    return PeakMemoryUsage() - m_max_bytes_used;
}

std::ostream& operator<<(std::ostream& ostream, MemorySystemQuery const& custom_memory_manager)
{
    ostream << custom_memory_manager.Statistics();
    return ostream;
}

std::string MemorySystemQuery::Statistics(std::string const& caller) const
{
    std::shared_lock lock(mutex);
    std::ostringstream oss;
    if (!caller.empty())
    {
        oss << "<Caller : " << caller << ">\n";
    }
    oss << "Current memory query:"
        << "\n Total allocated bytes  : " << m_total_allocated_bytes
        << "\n Max bytes used         : " << m_max_bytes_used
        << '\n';
    return oss.str();
}

uint64_t MemorySystemQuery::CurrentMemoryUsage()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(),
                         reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc),
                         sizeof(pmc));
    return static_cast<uint64_t>(pmc.WorkingSetSize + pmc.PrivateUsage);
#else
    return MemoryMonitor::CurrentUsage();
#endif
}

uint64_t MemorySystemQuery::PeakMemoryUsage()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(),
                         &pmc,
                         sizeof(pmc));
    return static_cast<uint64_t>(pmc.PeakWorkingSetSize);
#else
    return MemoryMonitor::PeakUsage();
#endif
}
