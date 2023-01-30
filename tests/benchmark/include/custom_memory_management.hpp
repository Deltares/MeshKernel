#pragma once

#if (DO_MANAGE_MEMORY)

#include <cstdlib>
#include <malloc.h>
#include <new>

#include "benchmark_memory_manager.hpp"

inline static size_t MemoryBlockSize(void* ptr)
{
#if defined(_WIN32)
    return _msize(ptr);
#else
    return malloc_usable_size(ptr);
#endif
}
/// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_free(void* ptr)
{
    CUSTOM_MEMORY_MANAGER.Unregister(MemoryBlockSize(ptr));
    std::free(ptr);
    ptr = nullptr;
}
#define free(ptr) custom_free(ptr)

#if defined(_WIN32)
/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_aligned_free(void* ptr)
{
    CUSTOM_MEMORY_MANAGER.Unregister(MemoryBlockSize(ptr));
    _aligned_free(ptr);
    ptr = nullptr;
}
#define _aligned_free(ptr) custom_aligned_free(ptr)
#endif

/// @brief Custom std::malloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return Pointer to the beginning of newly allocated memory
inline static void* custom_malloc(size_t size)
{
    if (void* ptr = std::malloc(size))
    {
        CUSTOM_MEMORY_MANAGER.Register(size);
        return ptr;
    }
    return nullptr;
}
#define malloc(size) custom_malloc(size)

/// @brief Custom std::calloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[num] Number of objects
/// @param[size] Number of bytes of uninitialized storage to be allocated per object
/// @return Pointer to the beginning of newly allocated memory
inline static void* custom_calloc(std::size_t num, size_t size)
{
    if (void* ptr = std::calloc(num, size))
    {
        CUSTOM_MEMORY_MANAGER.Register(size * num);
        return ptr;
    }
    return nullptr;
}
#define calloc(num, size) custom_calloc(num, size)

/// @brief Custom std::realloc wrapper which registers the re-allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory area to be reallocated
/// @param [new_size] New size of the array
/// @return Pointer to the beginning of newly allocated memory
inline static void* custom_realloc(void* ptr, std::size_t new_size)
{
    if (void* new_ptr = std::realloc(ptr, new_size))
    {
        bool const ptr_was_freed = (new_ptr != ptr);
        CUSTOM_MEMORY_MANAGER.Register(new_size, ptr_was_freed);
        return new_ptr;
    }
    return nullptr;
}
#define realloc(ptr, new_size) custom_realloc(ptr, new_size)

/// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
///        which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return The pointer to the beginning of newly allocated memory
inline static void* custom_aligned_malloc(size_t size, size_t alignment)
{
#if defined(_WIN32)
    // std::aligned_alloc is not implemented in VS
    void* ptr = _aligned_malloc(size, alignment);
#else
    void* ptr = std::aligned_alloc(alignment, size);
#endif
    if (ptr)
    {
        CUSTOM_MEMORY_MANAGER.Register(size);
        return ptr;
    }
    return nullptr;
}
#if defined(_WIN32)
#define _aligned_malloc(size, alignment) custom_aligned_malloc(size, alignment)
#else
#define aligned_alloc(alignment, size) custom_aligned_malloc(size, alignment)
#endif

#endif // DO_MANAGE_MEMORY
