#include <iostream>

#include <benchmark/benchmark.h>

#include "memory_management.hpp"

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
    ::benchmark::RegisterMemoryManager(&CUSTOM_MEMORY_MANAGER);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return EXIT_FAILURE;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    return EXIT_SUCCESS;
}