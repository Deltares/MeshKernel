#pragma once

#include "custom_memory_manager.hpp"
#include "platform.hpp"

// Redefinition of xalloc/free  et al. functions

#define malloc(size) CUSTOM_MEMORY_MANAGER.Malloc(size)

#define calloc(num, size) CUSTOM_MEMORY_MANAGER.Calloc(num, size)

#define realloc(ptr, new_size) CUSTOM_MEMORY_MANAGER.Realloc(ptr, new_size)

#define free(ptr) CUSTOM_MEMORY_MANAGER.Free(ptr)

// std::aligned_alloc is part of the c++17 standard, however it's not available in MSVC due
// to C11 restrictions. It provides _aligned_malloc instead, which should be freed by _aligned_free.
// On the other hand, GCC provides std::aligned_alloc. std::free is used to free the associated memory.
// In the defs below, the call to std::aligned_alloc is redirected to custom_aligned_alloc
// but under Linux, the order of parameters is reversed in order to have a unified function signature.
// free frees ptr allocated by aligned_alloc.
#if defined(WIN_MSVC)
#define _aligned_malloc(size, alignment) CUSTOM_MEMORY_MANAGER.AlignedAlloc(size, alignment)
#define _aligned_free(ptr) CUSTOM_MEMORY_MANAGER.AlignedFree(ptr)
#elif defined(LINUX_GNUC)
#define aligned_alloc(alignment, size) CUSTOM_MEMORY_MANAGER.AlignedAlloc(size, alignment)
#endif
