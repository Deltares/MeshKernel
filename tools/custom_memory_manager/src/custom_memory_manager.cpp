#include "custom_memory_manager.hpp"

#include <algorithm>
#include <sstream>

#include <malloc.h>

#include "platform.hpp"

CustomMemoryManager& CustomMemoryManager::Instance()
{
    static CustomMemoryManager instance;
    return instance;
}

void CustomMemoryManager::Start()
{
    ResetStatistics();
}

void CustomMemoryManager::Stop(Result* result)
{
    std::unique_lock lock(mutex);
    result->num_allocs = m_num_allocations;
    result->max_bytes_used = m_max_bytes_used;
    result->total_allocated_bytes = m_total_allocated_bytes;
    result->net_heap_growth = m_net_heap_growth;
}

void* CustomMemoryManager::Alloc(size_t size)
{
    std::unique_lock lock(mutex);
    if (void* ptr = std::malloc(size))
    {
        Register(size);
        return ptr;
    }
    return nullptr;
}

void* CustomMemoryManager::AlignedAlloc(size_t size, size_t alignment)
{
    std::unique_lock lock(mutex);
    void* ptr = nullptr;
#if defined(WIN_MSVC_BENCHMARK)
    // std::aligned_alloc is not implemented in VS
    ptr = _aligned_malloc(size, alignment);
#elif defined(LINUX_OR_APPLE_GNUC_BENCHMARK)
    ptr = std::aligned_alloc(alignment, size);
#endif
    if (ptr)
    {
        Register(size);
        return ptr;
    }
    return nullptr;
}

void CustomMemoryManager::Free(void* ptr)
{
    std::unique_lock lock(mutex);
    Unregister(MemoryBlockSize(ptr));
    std::free(ptr);
    ptr = nullptr;
}

void CustomMemoryManager::AlignedFree(void* ptr)
{
    std::unique_lock lock(mutex);
    Unregister(MemoryBlockSize(ptr));
#if defined(WIN_MSVC_BENCHMARK)
    _aligned_free(ptr);
#elif defined(LINUX_OR_APPLE_GNUC_BENCHMARK)
    std::free(ptr);
#endif
    ptr = nullptr;
}

void CustomMemoryManager::ResetStatistics()
{
    std::unique_lock lock(mutex);
    m_num_allocations = 0;
    m_num_deallocations = 0;
    m_total_allocated_bytes = 0;
    m_max_bytes_used = TombstoneValue;
    m_net_heap_growth = 0;
}

std::string CustomMemoryManager::Statistics(std::string const& caller) const
{
    std::shared_lock lock(mutex);
    std::ostringstream oss;
    if (!caller.empty())
    {
        oss << "<Caller : " << caller << ">\n";
    }
    oss << "Current memory manager statistics:"
        << "\n Number of allocations  : " << m_num_allocations
        << "\n Number of deallocations: " << m_num_deallocations
        << "\n Total allocated bytes  : " << m_total_allocated_bytes
        << "\n Max bytes used         : " << m_max_bytes_used
        << "\n Net heap growth        : " << m_net_heap_growth
        << '\n';
    return oss.str();
}

int64_t CustomMemoryManager::Allocations() const
{
    std::shared_lock lock(mutex);
    return m_num_allocations;
}

int64_t CustomMemoryManager::Deallocations() const
{
    std::shared_lock lock(mutex);
    return m_num_deallocations;
}

int64_t CustomMemoryManager::TotalAllocatedBytes() const
{
    std::shared_lock lock(mutex);
    return m_total_allocated_bytes;
}

int64_t CustomMemoryManager::MaxBytesUsed() const
{
    std::shared_lock lock(mutex);
    return m_max_bytes_used;
}

int64_t CustomMemoryManager::NetHeapGrowth() const
{
    std::shared_lock lock(mutex);
    return m_net_heap_growth;
}

bool CustomMemoryManager::HasLeaks() const
{
    return (Allocations() != Deallocations()) || (NetHeapGrowth() != 0);
}

void CustomMemoryManager::Register(int64_t size)
{
    m_num_allocations++;
    m_total_allocated_bytes += size;
    m_net_heap_growth += size;
    m_max_bytes_used = std::min(m_max_bytes_used, m_net_heap_growth);
}

void CustomMemoryManager::Unregister(int64_t size)
{
    m_num_deallocations++;
    m_net_heap_growth -= static_cast<int64_t>(size);
    m_max_bytes_used = std::max(m_max_bytes_used, m_net_heap_growth);
}

size_t CustomMemoryManager::MemoryBlockSize(void* ptr)
{
    if (ptr)
    {
#if defined(WIN_MSVC_BENCHMARK)
        return _msize(ptr);
#elif defined(LINUX_OR_APPLE_GNUC_BENCHMARK)
        return malloc_usable_size(ptr);
#endif
    }
    else
    {
        return 0;
    }
}

std::ostream& operator<<(std::ostream& ostream, CustomMemoryManager const& custom_memory_manager)
{
    ostream << custom_memory_manager.Statistics();
    return ostream;
}
