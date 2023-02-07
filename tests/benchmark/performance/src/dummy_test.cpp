#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include "memory_management.hpp"

static std::vector<double> VectorFilledWithRandomNumbers(size_t n, double lower_bound, double upper_bound)
{
    std::random_device random_device;
    std::mt19937 mersenne_engine{random_device()};
    std::normal_distribution<double> distribution(lower_bound, upper_bound);

    auto generator = [&distribution, &mersenne_engine]()
    {
        return distribution(mersenne_engine);
    };

    std::vector<double> vector(n);
    generate(vector.begin(), vector.end(), generator);
    return vector;
}

static double SumOfPowerOfVectorElements(std::vector<double> const& vector, double n)
{
    double result{0.0};
    std::for_each(vector.begin(), vector.end(), [&](auto const element)
                  { result += std::pow(element, n); });
    return result;
}

static void BM_Test1(benchmark::State& state)
{
    for (auto _ : state)
    {
        auto const vector = VectorFilledWithRandomNumbers(10000, 1, 10);
        SumOfPowerOfVectorElements(vector, 4);
    }
}
BENCHMARK(BM_Test1);

static void BM_Test2(benchmark::State& state)
{
    for (auto _ : state)
    {
        auto const vector = VectorFilledWithRandomNumbers(100000, 1, 100);
        SumOfPowerOfVectorElements(vector, 5.666);
    }
}
BENCHMARK(BM_Test2);

static void BM_Test2Cycles(benchmark::State& state)
{
    for (auto _ : state)
    {
        for (int i = 0; i < 10; ++i)
        {
            auto vector = VectorFilledWithRandomNumbers(100000, 1, 100);
            SumOfPowerOfVectorElements(vector, 5.666);
            vector.clear();
            vector.shrink_to_fit();
        }
    }
}
BENCHMARK(BM_Test2Cycles);
