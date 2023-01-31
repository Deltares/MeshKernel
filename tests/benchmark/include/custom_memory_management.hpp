#pragma once

#if (DO_MANAGE_MEMORY)

#include <cstdlib>
#include <cstring>
#include <malloc.h>
#include <new>

#include "benchmark_memory_manager.hpp"

inline static size_t MemoryBlockSize(void* ptr)
{
#if defined(_WIN32) && defined(_MSC_VER)
    return _msize(ptr);
#elif defined(__linux__) && define(__GNUC__)
    return malloc_usable_size(ptr);
#else
#error "Unsupported platform and/or compiler"
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

#if defined(_WIN32)
#if defined(_MSC_VER)
/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
inline static void custom_aligned_free(void* ptr)
{
    CUSTOM_MEMORY_MANAGER.Unregister(MemoryBlockSize(ptr));
    _aligned_free(ptr);
    ptr = nullptr;
}
#else
#error "Unsupported compiler"
#endif
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
    if (!ptr)
    {
        // allocate a new memory block of new_size bytes (realloc is equivalent to malloc)
        if (ptr = std::malloc(new_size))
        {
            CUSTOM_MEMORY_MANAGER.Register(new_size);
            return ptr;
        }
        // malloc failed
        return nullptr;
    }
    else
    {
        if (new_size == 0)
        {
            // std::free(ptr);
            // CUSTOM_MEMORY_MANAGER.Unregister(old_size);
            return nullptr;
        }
        else if (new_size <= MemoryBlockSize(ptr))
        {
            // memory block is sufficiently large, do nothing
            return ptr;
        }
        else
        {
            // first, try to expand without copying
#if defined(_WIN32) && defined(_MSC_VER)
            void* new_ptr = _expand(ptr, new_size);
#elif defined(__linux__) && define(__GNUC__)
            // TODO: what's the GCC equiv?
#error "Equivalent of _expand for GNU is needed"
#else
#error "Unsupported platform and/or compiler"
#endif
            if (new_ptr)
            {
                std::cout << "called expand, address: old = " << ptr << " new = " << new_ptr << '\n';
                CUSTOM_MEMORY_MANAGER.Register(new_size - MemoryBlockSize(ptr), false);
                return new_ptr;
            }
            else
            {
                // expansion failed, allocate a new memory block of new_size bytes
                if (new_ptr = std::malloc(new_size))
                {
                    std::cout << "called malloc\n";
                    // malloc succeeded, copy the old memory to the new one then free the old memory
                    size_t const old_size = MemoryBlockSize(ptr);
                    std::cout << "realloc size   : " << old_size << ' ' << new_size << '\n';
                    std::cout << "realloc address: " << ptr << ' ' << new_ptr << '\n';
                    std::memcpy(new_ptr, ptr, old_size);
                    std::free(ptr);
                    ptr = nullptr;
                    CUSTOM_MEMORY_MANAGER.Register(new_size);
                    CUSTOM_MEMORY_MANAGER.Unregister(old_size);
                    return new_ptr;
                }
                // malloc failed, a nullptr is returned, the old memory remains unchanged
                return nullptr;
            }
        }
    }
}

/// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
///        which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return The pointer to the beginning of newly allocated memory
inline static void* custom_aligned_malloc(size_t size, size_t alignment)
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

// Redefnition of xalloc/free functions (must be done after all custom functions are defined)
// realloc uses

#define free(ptr) custom_free(ptr)

#if defined(_WIN32)
#if defined(_MSC_VER)
#define _aligned_free(ptr) custom_aligned_free(ptr)
#else
#error "Unsupported compiler"
#endif
#endif

#define malloc(size) custom_malloc(size)

#define calloc(num, size) custom_calloc(num, size)

#define realloc(ptr, new_size) custom_realloc(ptr, new_size)

#if defined(_WIN32) && defined(_MSC_VER)
#define _aligned_malloc(size, alignment) custom_aligned_malloc(size, alignment)
#elif defined(__linux__) && define(__GNUC__) // this is really part of c++17
#define aligned_alloc(alignment, size) custom_aligned_malloc(size, alignment)
#else
#error "Unsupported platform and/or compiler"
#endif

#endif // DO_MANAGE_MEMORY
