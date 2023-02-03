#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <benchmark/benchmark.h>

#include "memory_management.hpp"

static void BM_memory1(benchmark::State& state)
{
    for (auto _ : state)
    {
        for (int i = 0; i < 10; ++i)
        {
            std::shared_ptr<std::vector<double>> ptr(std::make_shared<std::vector<double>>());
            ::benchmark::DoNotOptimize(ptr);
            for (size_t j = 0; j < 10; ++j)
            {
                ptr->push_back(static_cast<double>(i + j));
            }
        }
    }
}
BENCHMARK(BM_memory1);

static void BM_memory2(benchmark::State& state)
{
    for (auto _ : state)
    {
        for (int i = 0; i < 10; ++i)
        {
            std::vector<std::shared_ptr<size_t>> vector;
            ::benchmark::DoNotOptimize(vector);
            for (size_t j = 0; j < 100; ++j)
            {
                vector.push_back(std::make_shared<size_t>());
            }
            vector.clear();
            vector.shrink_to_fit();
        }
    }
}
BENCHMARK(BM_memory2);
