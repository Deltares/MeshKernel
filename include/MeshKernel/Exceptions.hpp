#pragma once
#include <stdexcept>

#include <MeshKernel/Constants.hpp>

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

    /// @brief Custom exception to describe an error caused by an invalid mesh
    class MeshGeometryError : public std::runtime_error
    {
    public:
        MeshGeometryError() : runtime_error(""), m_invalidIndex(sizetMissingValue), m_location(MeshLocations::None){};

        /// @brief Exception for mesh geometry errors accepting a string, the invalid mesh location index, and the location type
        explicit MeshGeometryError(const std::string& msg, size_t invalidIndex, MeshLocations location) : runtime_error(msg), m_invalidIndex(invalidIndex), m_location(location)
        {
        }
        /// @brief Exception for algorithmic errors accepting a char array, the invalid mesh location index, and the location type
        explicit MeshGeometryError(const char* msg, size_t invalidIndex, MeshLocations location) : runtime_error(msg), m_invalidIndex(invalidIndex), m_location(location) {}

        size_t m_invalidIndex;    ///< The invalid mesh location index
        MeshLocations m_location; ///< The location type
    };

} // namespace meshkernel
