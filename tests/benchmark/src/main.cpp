#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <new>
#include <vector>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

/// @brief Custom memory manager for registration via ::benchmark::RegisterMemoryManager
class BenchmarkMemoryManager final : public benchmark::MemoryManager
{
public:
    static BenchmarkMemoryManager& Instance()
    {
        static BenchmarkMemoryManager instance;
        return instance;
    }

    /// @brief Starts recording allocation information
    void Start() override
    {
        m_num_allocs = 0;
        m_num_deallocs = 0;
        m_total_allocated_bytes = 0;
        m_max_bytes_used = 0;
    }

    /// @brief Stop recording and fills out the given Result structure
    /// @param[result] Structure
    void Stop(Result* result) override
    {
        result->num_allocs = m_num_allocs;
        result->total_allocated_bytes = m_total_allocated_bytes;
        result->max_bytes_used = m_max_bytes_used;
    }

    /// @brief Registers an allocation
    /// @param[size] size of allocated memory
    void RegisterAllocation(size_t size)
    {
        m_num_allocs++;
        m_total_allocated_bytes += size;
        m_max_bytes_used = std::max(m_max_bytes_used, m_total_allocated_bytes);
    }

    /// @brief Registers an deallocation
    /// @param[ptr] size of allocated memory
    void RegisterDeallocation(void const* const ptr)
    {
        m_num_deallocs++;
        m_total_allocated_bytes -= sizeof(ptr); // this is incorrect!!!
    }

    size_t Allocations() const { return m_num_allocs; }

    size_t Deallocations() const { return m_num_deallocs; }

private:
    BenchmarkMemoryManager() = default;
    ~BenchmarkMemoryManager() = default;
    BenchmarkMemoryManager(BenchmarkMemoryManager const&) = delete;
    BenchmarkMemoryManager& operator=(BenchmarkMemoryManager const&) = delete;

    int64_t m_num_allocs = 0;            ///< The number of allocations made in total between Start and Stop
    int64_t m_num_deallocs = 0;          ///< The number of deallocations made in total between Start and Stop
    int64_t m_total_allocated_bytes = 0; ///< The total memory allocated, in bytes, between Start and Stop
    int64_t m_max_bytes_used = 0;        ///< The peak memory use between Start and Stop
};

#define MEMORY_MANAGER BenchmarkMemoryManager::Instance()

/// @brief Custom std::free wrapper which registers the deallocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
static void custom_free(void* ptr)
{
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    std::free(ptr);
}
#define free(ptr) custom_free(ptr)

/// @brief Custom std::malloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @return Pointer to the beginning of newly allocated memory
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

/// @brief Custom std::calloc wrapper which registers the allocation in a global BenchmarkMemoryManager object
/// @param[num] Number of objects
/// @param[size] Number of bytes of uninitialized storage to be allocated per object
/// @return Pointer to the beginning of newly allocated memory
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

/// @brief Custom std::realloc wrapper which registers the re-allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory area to be reallocated
/// @param [new_size] New size of the array
/// @return Pointer to the beginning of newly allocated memory
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

/// @brief Custom _aligned_malloc (WIN32) or std::aligned_alloc (LINUX) wrapper
///        which registers the allocation in a global BenchmarkMemoryManager object
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return The pointer to the beginning of newly allocated memory
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
/// @brief Custom _aligned_free wrapper (WIN32 only, free does the job under LINUX)
///        and registers the allocation in a global BenchmarkMemoryManager object
/// @param[ptr] Pointer to the memory block to deallocate
static void custom_aligned_free(void* ptr)
{
    MEMORY_MANAGER.RegisterDeallocation(ptr);
    _aligned_free(ptr);
}
#define _aligned_free(ptr) custom_aligned_free(ptr)
#endif

// Note on gloabl replacement of the new operator:
// Replacing the throwing single object allocation functions is sufficient to handle all allocations.
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_new

/// @brief Gloabl replacement (oveload) for void* operator new ( std::size_t size )
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

/// @brief Gloabl replacement (overload) for void* operator new ( std::size_t size, std::align_val_t alignment)
/// @param[size] Number of bytes of uninitialized storage to be allocated
/// @param[alignment] Specifies the alignment
/// @return If successful, a non-null pointer pointing to aligned memory, throws an allocation failure exception otherwise
void* operator new(std::size_t size, std::align_val_t alignment)
{
#if defined(WIN32)
    auto ptr = _aligned_malloc(size, static_cast<std::size_t>(alignment));
#else
    auto ptr = aligned_alloc(static_cast<std::size_t>(alignment), size);
#endif
    if (ptr)
    {
        return ptr;
    }
    throw std::bad_alloc{};
}

// Note on gloabl replacement of the delete operator:
// Rreplacing the throwing single object deallocation functions is sufficient to handle all deallocation overloads
// See section "Global replacements" in https://en.cppreference.com/w/cpp/memory/new/operator_delete

/// @brief Gloabl replacement (oveload) for void operator delete ( void* ptr ) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
void operator delete(void* ptr) noexcept { free(ptr); }

/// @brief Gloabl replacement (oveload) for void operator delete ( void* ptr ) noexcept
/// @param[ptr] Pointer to the memory block to deallocate
/// @param[alignment] Specifies the alignment (unused)
void operator delete(void* ptr, std::size_t size, std::align_val_t /*alignment*/) noexcept
{
#if defined(WIN32)
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

struct Pt
{
    int x = -1;
    int y = -1;
};

TEST(MemoryManager, new_operator)
{
    int* pt1 = new int[5];
    delete[] pt1;
    BenchmarkMemoryManager const& memory_manager = MEMORY_MANAGER;
    EXPECT_EQ(memory_manager.Allocations(), memory_manager.Deallocations());
}

int main(int argc, char** argv)
{

    int* pt1 = new int[5];
    delete[] pt1;
    std::vector<Pt> vec;
    BenchmarkMemoryManager const& memory_manager = MEMORY_MANAGER;
    std::cout << memory_manager.Allocations() << ' ' << memory_manager.Deallocations() << std::endl;

    ::testing::InitGoogleTest(&argc, argv);
    int test_ret = RUN_ALL_TESTS();

    if (!argv)
    {
        argc = 1;
        char arg0_default[] = "benchmark";
        char* args_default = arg0_default;
        argv = &args_default;
    }
    ::benchmark::Initialize(&argc, argv);

    ::benchmark::SetDefaultTimeUnit(::benchmark::kMillisecond);
    ::benchmark::RegisterMemoryManager(&MEMORY_MANAGER);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return EXIT_FAILURE;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return test_ret;
}