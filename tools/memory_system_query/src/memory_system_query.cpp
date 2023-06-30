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
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fstream>
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
    m_total_allocated_bytes = GetRAMSystemUsedByCurrentProcess();
    m_max_bytes_used = GetRAMPhysicalUsedByCurrentProcessPeak();
}

void MemorySystemQuery::Stop(Result* result)
{
    std::unique_lock lock(mutex);
    result->total_allocated_bytes = GetRAMSystemUsedByCurrentProcess() - m_total_allocated_bytes;
    result->max_bytes_used = GetRAMPhysicalUsedByCurrentProcessPeak() - m_max_bytes_used;
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

//
// #ifdef __linux__
// namespace
//{
//    int ParseLine(char* line)
//    {
//        const auto i = std::strlen(line);
//
//        while (*line < '0' || *line > '9')
//        {
//            line++;
//        }
//
//        line[i - 3] = '\0';
//        return std::atoi(line);
//    }
//
//    std::string const proc_self_status_path("/proc/self/status");
//
//    int64_t KiloBytesToBytes(int64_t kilo_bytes)
//    {
//        return kilo_bytes * 1024;
//    }
//
//    int64_t Physical()
//    {
//        int64_t kilo_bytes = 0;
//        std::ifstream file(proc_self_status_path, std::ifstream::in);
//        std::string line;
//        while (std::getline(file, line))
//        {
//            if (line.compare(1, 6, "VmRSS:") == 0)
//            {
//                kilo_bytes += ParseLine(line.data());
//            }
//        }
//        file.close();
//        return KiloBytesToBytes(kilo_bytes);
//    }
//
//    int64_t PhysicalPeak()
//    {
//        int64_t kilo_bytes = 0;
//        std::ifstream file(proc_self_status_path, std::ifstream::in);
//        std::string line;
//        while (std::getline(file, line))
//        {
//            if (line.compare(1, 6, "VmHWM:") == 0)
//            {
//                kilo_bytes += ParseLine(line.data());
//            }
//        }
//        file.close();
//        return KiloBytesToBytes(kilo_bytes);
//    }
//
//    int64_t Virtual()
//    {
//        int64_t kilo_bytes = 0;
//        std::ifstream file(proc_self_status_path, std::ifstream::in);
//        std::string line;
//        while (std::getline(file, line))
//        {
//            if (line.compare(1, 7, "VmSize:") == 0)
//            {
//                kilo_bytes = ParseLine(line.data());
//                break;
//            }
//        }
//        file.close();
//        return KiloBytesToBytes(kilo_bytes);
//    }
//} // namespace
// #endif
//
// int64_t MemorySystemQuery::GetRAMSystemUsedByCurrentProcess()
//{
// #ifdef _WIN32
//    PROCESS_MEMORY_COUNTERS_EX pmc;
//    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
//    return static_cast<int64_t>(pmc.WorkingSetSize);
// #else
//    return GetRAMPhysicalUsedByCurrentProcess() + GetRAMVirtualUsedByCurrentProcess();
// #endif
//}
//
// int64_t MemorySystemQuery::GetRAMPhysicalUsedByCurrentProcess()
//{
// #ifdef _WIN32
//    PROCESS_MEMORY_COUNTERS_EX pmc;
//    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
//    return static_cast<int64_t>(pmc.WorkingSetSize);
// #else
//    return Physical();
// #endif
//}
//
// int64_t MemorySystemQuery::GetRAMPhysicalUsedByCurrentProcessPeak()
//{
// #if defined(_WIN32)
//    PROCESS_MEMORY_COUNTERS pmc;
//    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
//    return static_cast<int64_t>(pmc.PeakWorkingSetSize);
// #else
//    return PhysicalPeak();
// #endif
//}
//
// int64_t MemorySystemQuery::GetRAMVirtualUsedByCurrentProcess()
//{
// #ifdef _WIN32
//    PROCESS_MEMORY_COUNTERS_EX pmc;
//    GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PPROCESS_MEMORY_COUNTERS>(&pmc), sizeof(pmc));
//    return pmc.PrivateUsage;
// #else
//    return Virtual();
// #endif
//}
