#pragma once

#if (DO_MANAGE_MEMORY)

#include <cstdlib>

#include "benchmark_memory_manager.hpp"

/// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_free(void* ptr)
{
    // printf("custom_free\n");
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    std::free(ptr);
}
#define free(ptr) custom_free(ptr)

/// @brief Custom std::malloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return Pointer to the beginning of newly allocated memory
inline static void* custom_malloc(size_t size)
{
    // printf("custom_malloc\n");
    if (size == 0)
    {
        return nullptr;
    }
    if (void* ptr = std::malloc(size))
    {
        MEMORY_MANAGER.RegisterAllocation(ptr, size);
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
    // printf("custom_calloc\n");
    if (size == 0)
    {
        return nullptr;
    }
    if (void* ptr = std::calloc(num, size))
    {
        MEMORY_MANAGER.RegisterAllocation(ptr, size * num);
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
    // printf("custom_realloc\n");
    if (new_size == 0)
    {
        free(ptr); // should this be this function's responsibility? Prob not.
        return nullptr;
    }
    if (void* new_ptr = std::realloc(ptr, new_size))
    {
        MEMORY_MANAGER.RegisterAllocation(new_ptr, new_size);
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
    // printf("custom_aligned_malloc\n");
#if defined(_WIN32)
    // std::aligned_alloc is not implemented in VS
    void* ptr = _aligned_malloc(size, alignment);
#else
    void* ptr = std::aligned_alloc(alignment, size);
#endif
    if (ptr)
    {
        MEMORY_MANAGER.RegisterAllocation(ptr, size);
        return ptr;
    }
    return nullptr;
}
#if defined(_WIN32)
#define _aligned_malloc(size, alignment) custom_aligned_malloc(size, alignment)
#else
#define aligned_alloc(alignment, size) custom_aligned_malloc(size, alignment)
#endif

#if defined(_WIN32)
/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_aligned_free(void* ptr)
{
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    _aligned_free(ptr);
}
#define _aligned_free(ptr) custom_aligned_free(ptr)
#endif

#endif // DO_MANAGE_MEMORY
