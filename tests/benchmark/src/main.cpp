#include <iostream>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

// #define DO_MANAGE_MEMORY true
#include "benchmark_memory_manager.hpp"
// #include "custom_memory_management.cpp"
#include "custom_memory_management.hpp"

int main(int argc, char** argv)
{
    // Unit tests
    ::testing::InitGoogleTest(&argc, argv);
    if (RUN_ALL_TESTS() != 0)
    {
        return EXIT_FAILURE;
    }

    // Benchmarks
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