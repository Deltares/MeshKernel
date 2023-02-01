#pragma once

#include <cstdlib>
#include <cstring>
#include <malloc.h>

#include "benchmark_memory_manager.hpp"

/// @brief Returns the size, in bytes, of a memory block allocated in the heap
/// @param[ptr] Pointer to the memory block
inline static size_t MemoryBlockSize(void* ptr)
{
    if (ptr)
    {
#if defined(_WIN32) && defined(_MSC_VER)
        return _msize(ptr);
#elif defined(__linux__) && define(__GNUC__)
        return malloc_usable_size(ptr);
#else
#error "Unsupported platform and/or compiler"
#endif
    }
    else
    {
        return 0;
    }
}

/// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_free(void* ptr)
{
    CUSTOM_MEMORY_MANAGER.Unregister(MemoryBlockSize(ptr));
    std::free(ptr);
    ptr = nullptr;
}

/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_aligned_free(void* ptr)
{
    CUSTOM_MEMORY_MANAGER.Unregister(MemoryBlockSize(ptr));
#if defined(_WIN32) && defined(_MSC_VER)
    _aligned_free(ptr);
#elif defined(__linux__) && define(__GNUC__)
    std::free(ptr);
#else
#error "Unsupported platform and/or compiler"
#endif
    ptr = nullptr;
}

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

/// @brief Custom std::realloc wrapper which registers the re-allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory area to be reallocated
/// @param [new_size] New size of the array
/// @return Pointer to the beginning of newly allocated memory
inline static void* custom_realloc(void* ptr, std::size_t new_size)
{
    // store the old size prior toreallocation
    size_t const old_size = MemoryBlockSize(ptr);
    if (void* new_ptr = std::realloc(ptr, new_size))
    {
        if (new_ptr != ptr)
        {
            // the address has changed because new ptr was malloced
            // register new allocation due malloc of new ptr
            CUSTOM_MEMORY_MANAGER.Register(new_size);
            // register deallocation due to free of old ptr
            CUSTOM_MEMORY_MANAGER.Unregister(old_size);
        }
        else
        {
            // the address did not change, the memory block was either expanded or shrunk
            // register only the size difference without incrementing the number of allocations
            CUSTOM_MEMORY_MANAGER.Register(new_size - old_size, false);
        }
        return new_ptr;
    }
    return nullptr;
}

/// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
///        which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return The pointer to the beginning of newly allocated memory
inline static void* custom_aligned_alloc(size_t size, size_t alignment)
{
#if defined(_WIN32) && defined(_MSC_VER)
    // std::aligned_alloc is not implemented in VS
    void* ptr = _aligned_malloc(size, alignment);
#elif defined(__linux__) && define(__GNUC__)
    void* ptr = std::aligned_alloc(alignment, size);
#else
#error "Unsupported platform and/or compiler"
#endif
    if (ptr)
    {
        CUSTOM_MEMORY_MANAGER.Register(size);
        return ptr;
    }
    return nullptr;
}

// Redefinition of xalloc/free  et al. functions

#define malloc(size) custom_malloc(size)

#define calloc(num, size) custom_calloc(num, size)

#define realloc(ptr, new_size) custom_realloc(ptr, new_size)

#define free(ptr) custom_free(ptr)

// std::aligned_alloc is part of c++17 standard, however it's not available in MSVC due
// to C11 restrictions. It provides _aligned_malloc instead, which should be freed by _aligned_free.
// On the other hand, GCC provides std::aligned_alloc. std::free is used to free the associated memory.
// In the defs below, the call to std::aligned_alloc is redirected to custom_aligned_alloc
// but under Linux, the order of parameters is reversed in order to have a unified function signature.
#if defined(_WIN32) && defined(_MSC_VER)
#define _aligned_malloc(size, alignment) custom_aligned_alloc(size, alignment)
#define _aligned_free(ptr) custom_aligned_free(ptr)
#elif defined(__linux__) && define(__GNUC__)
#define aligned_alloc(alignment, size) custom_aligned_alloc(size, alignment)
#error "Unsupported platform and/or compiler"
#endif
