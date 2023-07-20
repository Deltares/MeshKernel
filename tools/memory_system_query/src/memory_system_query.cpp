#include "memory_system_query.hpp"

#include <sstream>
#include <string>

MemorySystemQuery& MemorySystemQuery::Instance()
{
    static MemorySystemQuery instance;
    return instance;
}

void MemorySystemQuery::Start()
{
    std::unique_lock lock(m_mutex);
    m_total_allocated_bytes = CurrentMemoryUsage();
    m_max_bytes_used = PeakMemoryUsage();
    /*std::cout << "start: "
              << "phys: " << m_memory_monitor->Usage(MemoryType::Physical)
              << " virt: " << m_memory_monitor->Usage(MemoryType::Virtual) << std::endl;*/
}

void MemorySystemQuery::Stop(Result& result)
{
    std::unique_lock lock(m_mutex);
    result.total_allocated_bytes = CurrentMemoryUsage() - m_total_allocated_bytes;
    result.max_bytes_used = PeakMemoryUsage() - m_max_bytes_used;
    /*std::cout << "end  : "
              << "phys: " << m_memory_monitor->Usage(MemoryType::Physical)
              << " virt: " << m_memory_monitor->Usage(MemoryType::Virtual) << std::endl;
    std::cout << "diff : " << result.total_allocated_bytes << std::endl;*/
}

int64_t MemorySystemQuery::TotalAllocatedBytes() const
{
    std::shared_lock lock(m_mutex);
    return CurrentMemoryUsage() - m_total_allocated_bytes;
}

int64_t MemorySystemQuery::MaxBytesUsed() const
{
    std::shared_lock lock(m_mutex);
    return PeakMemoryUsage() - m_max_bytes_used;
}

std::string MemorySystemQuery::Statistics(std::string const& caller) const
{
    std::shared_lock lock(m_mutex);
    std::ostringstream oss;
    if (!caller.empty())
    {
        oss << "<Caller : " << caller << ">\n";
    }
    oss << "Current memory query:"
        << "\n Total allocated bytes : " << TotalAllocatedBytes()
        << "\n Max bytes used        : " << MaxBytesUsed()
        << '\n';
    return oss.str();
}

uint64_t MemorySystemQuery::CurrentMemoryUsage() const
{
    return m_memory_monitor->CurrentUsage();
}

uint64_t MemorySystemQuery::PeakMemoryUsage() const
{
    return m_memory_monitor->PeakUsage();
}

std::ostream& operator<<(std::ostream& ostream,
                         MemorySystemQuery const& custom_memory_manager)
{
    ostream << custom_memory_manager.Statistics();
    return ostream;
}