#include <new>

#include <benchmark/benchmark.h>

class BenchmarkMemoryManager : public benchmark::MemoryManager
{
public:
    int64_t num_allocs;
    int64_t max_bytes_used;

    void Start() override
    {
        num_allocs = 0;
        max_bytes_used = 0;
    }

    void Stop(Result* result) override
    {
        result->num_allocs = num_allocs;
        result->max_bytes_used = max_bytes_used;
    }
};

static BenchmarkMemoryManager memory_manager;

// overload new to allow increasing the number of allocations and num of max bytes used (peak mem)
void* operator new(size_t size)
{
    // std::cout << " >> custom new" << std::endl;
    if (size == 0)
    {
        // avoid std::malloc(0) which may return nullptr on success
        ++size;
    }
    if (void* ptr = std::malloc(size))
    {
        memory_manager.num_allocs++;
        memory_manager.max_bytes_used += size;
        return ptr;
    }
    throw std::bad_alloc{};
}

int main(int argc, char** argv)
{
    ::benchmark::SetDefaultTimeUnit(::benchmark::kMillisecond);
    ::benchmark::RegisterMemoryManager(&memory_manager);
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::RegisterMemoryManager(nullptr);
}