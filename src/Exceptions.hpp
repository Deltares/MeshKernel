#pragma once

#include <stdexcept>

namespace meshkernel
{
    class AlgorithmError : public std::runtime_error
    {
    public:
        explicit AlgorithmError(const std::string& msg) : runtime_error(msg) {}
        explicit AlgorithmError(const char* msg) : runtime_error(msg) {}
    };
} // namespace meshkernel
