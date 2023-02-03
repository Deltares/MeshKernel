#include "custom_memory_manager.hpp"

#include <algorithm>
#include <malloc.h>
#include <sstream>

#include "platform.hpp"

CustomMemoryManager& CustomMemoryManager::Instance()
{
    static CustomMemoryManager instance;
    return instance;
}

void CustomMemoryManager::Start() { ResetStatistics(); }

void CustomMemoryManager::Stop(Result* result)
{
    result->num_allocs = m_num_allocations;
    result->max_bytes_used = m_max_bytes_used;
    result->total_allocated_bytes = m_total_allocated_bytes;
    result->net_heap_growth = m_net_heap_growth;
}

void CustomMemoryManager::Free(void* ptr)
{
    Unregister(MemoryBlockSize(ptr));
    std::free(ptr);
    ptr = nullptr;
}

void CustomMemoryManager::AlignedFree(void* ptr)
{
    Unregister(MemoryBlockSize(ptr));
#if defined(WIN_MSVC)
    _aligned_free(ptr);
#elif defined(LINUX_GNUC)
    std::free(ptr);
#endif
    ptr = nullptr;
}

void* CustomMemoryManager::Malloc(size_t size)
{
    if (void* ptr = std::malloc(size))
    {
        Register(size);
        return ptr;
    }
    return nullptr;
}

void* CustomMemoryManager::Calloc(std::size_t num, size_t size)
{
    if (void* ptr = std::calloc(num, size))
    {
        Register(size * num);
        return ptr;
    }
    return nullptr;
}

void* CustomMemoryManager::Realloc(void* ptr, std::size_t new_size)
{
    // store the old size prior toreallocation
    size_t const old_size = MemoryBlockSize(ptr);
    if (void* new_ptr = std::realloc(ptr, new_size))
    {
        if (new_ptr != ptr)
        {
            // the address has changed because new ptr was malloced
            // register new allocation due malloc of new ptr
            Register(new_size);
            // register deallocation due to free of old ptr
            Unregister(old_size);
        }
        else
        {
            // the address did not change, the memory block was either expanded or shrunk
            // register only the size difference without incrementing the number of allocations
            Register(new_size - old_size, false);
        }
        return new_ptr;
    }
    return nullptr;
}

void* CustomMemoryManager::AlignedAlloc(size_t size, size_t alignment)
{
    void* ptr = nullptr;
#if defined(WIN_MSVC)
    // std::aligned_alloc is not implemented in VS
    ptr = _aligned_malloc(size, alignment);
#elif defined(LINUX_GNUC)
    ptr = std::aligned_alloc(alignment, size);
#endif
    if (ptr)
    {
        Register(size);
        return ptr;
    }
    return nullptr;
}

void CustomMemoryManager::ResetStatistics()
{
    m_num_allocations = 0;
    m_num_deallocations = 0;
    m_total_allocated_bytes = 0;
    m_max_bytes_used = TombstoneValue;
    m_net_heap_growth = 0;
}

std::string CustomMemoryManager::Statistics(std::string const& caller) const
{
    std::ostringstream oss;
    if (!caller.empty())
    {
        oss << caller << '\n';
    }
    oss << *this;
    return oss.str();
}

int64_t CustomMemoryManager::Allocations() const { return m_num_allocations; }

int64_t CustomMemoryManager::Deallocations() const { return m_num_deallocations; }

int64_t CustomMemoryManager::TotalAllocatedBytes() const { return m_total_allocated_bytes; }

int64_t CustomMemoryManager::MaxBytesUsed() const { return m_max_bytes_used; }

int64_t CustomMemoryManager::NetHeapGrowth() const { return m_net_heap_growth; }

bool CustomMemoryManager::HasLeaks() const
{
    return ((m_num_allocations != m_num_deallocations) || (m_net_heap_growth != 0));
}

void CustomMemoryManager::Register(int64_t size, bool incrrement_num_allocations)
{
    if (incrrement_num_allocations)
    {
        m_num_allocations++;
    }
    m_total_allocated_bytes += size;
    m_net_heap_growth += size;
    m_max_bytes_used = std::min(m_max_bytes_used, m_net_heap_growth);
}

void CustomMemoryManager::Unregister(int64_t size)
{
    m_num_deallocations++;
    m_net_heap_growth -= size;
    m_max_bytes_used = std::max(m_max_bytes_used, m_net_heap_growth);
}

size_t CustomMemoryManager::MemoryBlockSize(void* ptr)
{
    if (ptr)
    {
#if defined(WIN_MSVC)
        return _msize(ptr);
#elif defined(LINUX_GNUC)
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
    ostream << "Current memory manager statistics:"
            << "\nNumber of allocations  : " << custom_memory_manager.m_num_allocations
            << "\nNumber of deallocations: " << custom_memory_manager.m_num_deallocations
            << "\nTotal allocated bytes  : " << custom_memory_manager.m_total_allocated_bytes
            << "\nMax bytes used         : " << custom_memory_manager.m_max_bytes_used
            << "\nNet heap growth        : " << custom_memory_manager.m_net_heap_growth
            << '\n';
    return ostream;
}
