#include <iostream>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

#include "benchmark_memory_manager.hpp"

#define DO_MANAGE_MEMORY true
#include "custom_memory_management.cpp"

int main(int argc, char** argv)
{
    // Unit tests
    ::testing::InitGoogleTest(&argc, argv);
    int test_ret = RUN_ALL_TESTS();

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
    ::benchmark::RegisterMemoryManager(&MEMORY_MANAGER);
    if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return EXIT_FAILURE;
    }
    ::benchmark::RunSpecifiedBenchmarks();
    ::benchmark::Shutdown();

    // is it possible to do an atexit test?
    std::cout << "MAIN TRACKMEM: "
              << MEMORY_MANAGER.Allocations() << ' '
              << MEMORY_MANAGER.Deallocations() << '\n';

    // return test_ret;

    return 0;
}