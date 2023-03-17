#include <benchmark/benchmark.h>

#ifdef ENABLE_MEM_REPORT
#include "custom_memory_manager.hpp"
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
#ifdef ENABLE_MEM_REPORT
    ::benchmark::RegisterMemoryManager(&CUSTOM_MEMORY_MANAGER);
#endif
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return EXIT_FAILURE;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return EXIT_SUCCESS;
}