#include "memory_system_query.hpp"

#include <map>
#include <sstream>
#include <string>
#include <utility>

// clang-format off
#if defined(_WIN32)
//
#include <windows.h>
#include <psapi.h>
//
#elif defined(__linux__)
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#else
//
#error "Unsupported platform"
//
#endif
// clang-format on

MemorySystemQuery& MemorySystemQuery::Instance()
{
    static MemorySystemQuery instance;
    return instance;
}

void MemorySystemQuery::Start()
{
    std::unique_lock lock(mutex);
    m_max_bytes_used = GetRAMPhysicalUsedByCurrentProcessPeak();
    m_total_allocated_bytes = GetRAMSystemUsedByCurrentProcess();
}

void MemorySystemQuery::Stop(Result* result)
{
    std::unique_lock lock(mutex);
    result->max_bytes_used = GetRAMPhysicalUsedByCurrentProcessPeak() - m_max_bytes_used;
    result->total_allocated_bytes = GetRAMSystemUsedByCurrentProcess() - m_total_allocated_bytes;
}

int64_t MemorySystemQuery::TotalAllocatedBytes() const
{
    std::shared_lock lock(mutex);
    return GetRAMSystemUsedByCurrentProcess() - m_total_allocated_bytes;
}

int64_t MemorySystemQuery::MaxBytesUsed() const
{
    std::shared_lock lock(mutex);
    return GetRAMPhysicalUsedByCurrentProcessPeak() - m_max_bytes_used;
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

int64_t MemorySystemQuery::GetRAMSystemUsedByCurrentProcess()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
    return static_cast<int64_t>(pmc.WorkingSetSize);
#else
    return GetRAMPhysicalUsedByCurrentProcess() + GetRAMVirtualUsedByCurrentProcess();
#endif
}

#ifdef __linux__

static int64_t constexpr Kilobytes2Bytes{1024};

static int ParseLine(char* line)
{
    const auto i = strlen(line);

    while (*line < '0' || *line > '9')
    {
        line++;
    }

    line[i - 3] = '\0';
    return atoi(line);
}

#endif

int64_t MemorySystemQuery::GetRAMPhysicalUsedByCurrentProcess()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
    return static_cast<int64_t>(pmc.WorkingSetSize);
#else
    int64_t result = 0;
    FILE* file = fopen("/proc/self/status", "r");
    constexpr int BufferSize{128};
    char line[BufferSize];
    while (fgets(line, BufferSize, file) != nullptr)
    {
        if (strncmp(line, "VmRSS:", 6) == 0)
        {
            result += ParseLine(line) * Kilobytes2Bytes;
        }
    }
    fclose(file);
    return result;
#endif
}

int64_t MemorySystemQuery::GetRAMPhysicalUsedByCurrentProcessPeak()
{
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
    return static_cast<int64_t>(pmc.PeakWorkingSetSize);
#else
    int64_t result = 0;
    FILE* file = fopen("/proc/self/status", "r");
    constexpr int BufferSize{128};
    char line[BufferSize];
    while (fgets(line, BufferSize, file) != nullptr)
    {
        if (strncmp(line, "VmHWM:", 6) == 0)
        {
            result += ParseLine(line) * Kilobytes2Bytes;
        }
    }
    fclose(file);
    return result;
#endif
}

int64_t MemorySystemQuery::GetRAMVirtualUsedByCurrentProcess()
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
    return pmc.PrivateUsage;
#else
    int64_t result = 0;
    constexpr int BufferSize{128};
    FILE* file = fopen("/proc/self/status", "r");
    char line[BufferSize];
    while (fgets(line, BufferSize, file) != nullptr)
    {
        if (strncmp(line, "VmSize:", 7) == 0)
        {
            result = ParseLine(line) * Kilobytes2Bytes;
            break;
        }
    }
    fclose(file);
    return result;
#endif
}
