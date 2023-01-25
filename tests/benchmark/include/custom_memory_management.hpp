#pragma once

#include "benchmark_memory_manager.hpp"

#if DO_MANAGE_MEMORY

/// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
static void custom_free(void* ptr);

/// @brief Custom std::malloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return Pointer to the beginning of newly allocated memory
static void* custom_malloc(size_t size);

/// @brief Custom std::calloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[num] Number of objects
/// @param[size] Number of bytes of uninitialized storage to be allocated per object
/// @return Pointer to the beginning of newly allocated memory
static void* custom_calloc(std::size_t num, size_t size);

/// @brief Custom std::realloc wrapper which registers the re-allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory area to be reallocated
/// @param [new_size] New size of the array
/// @return Pointer to the beginning of newly allocated memory
static void* custom_realloc(void* ptr, std::size_t new_size);

/// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
///        which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return The pointer to the beginning of newly allocated memory
static void* custom_aligned_malloc(size_t size, size_t alignment);

#if defined(WIN32)
/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
static void custom_aligned_free(void* ptr);
#endif

// Note on gloabl replacement of the new operator:
// Replacing the throwing single object allocation functions is sufficient to handle all allocations.
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_new

/// @brief Gloabl replacement (oveload) for void* operator new ( std::size_t size )
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return If successful, a non-null pointer, throws an allocation failure exception otherwise
void* operator new(size_t size);

/// @brief Gloabl replacement (overload) for void* operator new ( std::size_t size, std::align_val_t alignment)
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return If successful, a non-null pointer pointing to aligned memory, throws an allocation failure exception otherwise
void* operator new(std::size_t size, std::align_val_t alignment);

// Note on gloabl replacement of the delete operator:
// Rreplacing the throwing single object deallocation functions is sufficient to handle all deallocation overloads
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_delete

/// @brief Gloabl replacement (oveload) for void operator delete ( void* ptr ) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
void operator delete(void* ptr) noexcept;

/// @brief Gloabl replacement (oveload) for void operator delete ( void* ptr ) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
/// @param[alignment] Specifies the alignment (unused)
void operator delete(void* ptr, std::size_t size, std::align_val_t /*alignment*/) noexcept;

#endif // DO_MANAGE_MEMORY