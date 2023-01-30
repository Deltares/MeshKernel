#if (DO_MANAGE_MEMORY)

#include "custom_memory_management.hpp"

// Note on gloabl replacement for the new operator:
// Replacing the throwing single object allocation functions is sufficient to handle all allocations.
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_new

/// @brief Gloabl replacement for void* operator new ( std::size_t size )
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return If successful, a non-null pointer, throws an allocation failure exception otherwise
void* operator new(size_t size)
{
    if (void* ptr = malloc(size))
    {
        return ptr;
    }
    throw std::bad_alloc{};
}

/// @brief Gloabl replacement for void* operator new ( std::size_t size, std::align_val_t alignment)
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return If successful, a non-null pointer pointing to aligned memory, throws an allocation failure exception otherwise
void* operator new(std::size_t size, std::align_val_t alignment)
{
#if defined(_WIN32)
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

// Note on gloabl replacement for the delete operator:
// Rreplacing the throwing single object deallocation functions is sufficient to handle all deallocation
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_delete

/// @brief Gloabl replacement for void operator delete ( void* ptr ) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
void operator delete(void* ptr) noexcept { free(ptr); }

/// @brief Gloabl replacement for void operator delete(void* ptr, size_t size) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
void operator delete(void* ptr, size_t /*size*/) noexcept { free(ptr); }

/// @brief Gloabl replacement for void operator delete(void* ptr, std::align_val_t alignment) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
/// @param[alignment] Specifies the alignment (unused)
void operator delete(void* ptr, std::align_val_t /*alignment*/) noexcept
{
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

/// @brief Gloabl replacement for void operator delete(void* ptr, size_t size, std::align_val_t alignment) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
/// @param[alignment] Specifies the alignment (unused)
void operator delete(void* ptr, size_t /*size*/, std::align_val_t /*alignment*/) noexcept
{
#if defined(_WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

#endif // DO_MANAGE_MEMORY
