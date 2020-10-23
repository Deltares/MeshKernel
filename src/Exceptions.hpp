#pragma once

#include <stdexcept>

namespace meshkernel
{
    class AlgorithmError : public std::runtime_error
    {
    public:
        AlgorithmError(const std::string& msg) : runtime_error(msg) {}
        AlgorithmError(const char* msg) : runtime_error(msg) {}
    };
} // namespace meshkernel
