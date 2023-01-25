#pragma once

#include <cstdlib>
#include <new>

#include "custom_memory_management.hpp"

#if DO_MANAGE_MEMORY

static void custom_free(void* ptr)
{
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    std::free(ptr);
}
#define free(ptr) custom_free(ptr)

static void* custom_malloc(size_t size)
{
    if (size == 0)
    {
        return nullptr;
    }
    if (void* ptr = std::malloc(size))
    {
        MEMORY_MANAGER.RegisterAllocation(size);
        return ptr;
    }
    return nullptr;
}
#define malloc(size) custom_malloc(size)

static void* custom_calloc(std::size_t num, size_t size)
{
    if (size == 0)
    {
        return nullptr;
    }
    if (void* ptr = std::calloc(num, size))
    {
        MEMORY_MANAGER.RegisterAllocation(size * num);
        return ptr;
    }
    return nullptr;
}
#define calloc(num, size) custom_calloc(num, size)

static void* custom_realloc(void* ptr, std::size_t new_size)
{
    if (new_size == 0)
    {
        free(ptr); // should this be this function's responsibility? Prob not.
        return nullptr;
    }
    if (void* new_ptr = std::realloc(ptr, new_size))
    {
        MEMORY_MANAGER.RegisterAllocation(new_size); // Do we have to increase the number of allocations? Maybe...
        return new_ptr;
    }
    return nullptr;
}
#define realloc(ptr, new_size) custom_realloc(ptr, new_size)

static void* custom_aligned_malloc(size_t size, size_t alignment)
{
#if defined(WIN32)
    // std::aligned_alloc is not implemented in VS
    void* ptr = _aligned_malloc(size, alignment);
#else
    void* ptr = std::aligned_alloc(alignment, size);
#endif
    if (ptr)
    {
        MEMORY_MANAGER.RegisterAllocation(size);
        return ptr;
    }
    return nullptr;
}
#if defined(WIN32)
#define _aligned_malloc(size, alignment) custom_aligned_malloc(size, alignment)
#else
#define aligned_malloc(alignment, size) custom_aligned_malloc(size, alignment)
#endif

#if defined(WIN32)
static void custom_aligned_free(void* ptr)
{
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    _aligned_free(ptr);
}
#define _aligned_free(ptr) custom_aligned_free(ptr)
#endif

void* operator new(size_t size)
{
    if (void* ptr = malloc(size))
    {
        return ptr;
    }
    throw std::bad_alloc{};
}

void* operator new(std::size_t size, std::align_val_t alignment)
{
#if defined(WIN32)
    void* ptr = _aligned_malloc(size, static_cast<std::size_t>(alignment));
#else
    void* ptr = aligned_alloc(static_cast<std::size_t>(alignment), size);
#endif
    if (ptr)
    {
        return ptr;
    }
    throw std::bad_alloc{};
}

void operator delete(void* ptr) noexcept { free(ptr); }

void operator delete(void* ptr, std::size_t size, std::align_val_t /*alignment*/) noexcept
{
#if defined(WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

#endif // DO_MANAGE_MEMORY