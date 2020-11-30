#pragma once

#include <stdexcept>

namespace meshkernel
{
    /// @brief Custom exception to describe an error caused by an algorithm
    class AlgorithmError : public std::runtime_error
    {
    public:
        /// @brief Exception for algorithmic errors accepting a string
        explicit AlgorithmError(const std::string& msg) : runtime_error(msg) {}
        /// @brief Exception for algorithmic errors accepting a char array
        explicit AlgorithmError(const char* msg) : runtime_error(msg) {}
    };
} // namespace meshkernel
