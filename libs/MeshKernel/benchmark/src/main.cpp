#include <benchmark/benchmark.h>

#ifdef ENABLE_BENCHMARKING_MEM_REPORT
#if defined(MEM_COLLECTION_METHOD_COUNT_BYTES)
#include "custom_memory_manager.hpp"
#elif defined(MEM_COLLECTION_METHOD_QUERY_SYSTEM)
#include "memory_system_query.hpp"
#endif
#endif

int main(int argc, char** argv)
{
    if (!argv)
    {
        argc = 1;
        char arg0_default[] = "MeshKernelBenchmark";
        char* args_default = arg0_default;
        argv = &args_default;
    }
    ::benchmark::Initialize(&argc, argv);
    ::benchmark::SetDefaultTimeUnit(::benchmark::kMillisecond);
#ifdef ENABLE_BENCHMARKING_MEM_REPORT
#if defined(MEM_COLLECTION_METHOD_COUNT_BYTES)
    ::benchmark::RegisterMemoryManager(&CUSTOM_MEMORY_MANAGER);
#elif defined(MEM_COLLECTION_METHOD_QUERY_SYSTEM)
    ::benchmark::RegisterMemoryManager(&MEMORY_SYSTEM_QUERY);
#endif
#endif
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return EXIT_FAILURE;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return EXIT_SUCCESS;
}