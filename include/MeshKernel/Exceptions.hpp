#pragma once

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Mesh.hpp>

namespace meshkernel
{
    /// @brief Custom exception to describe an error caused by an algorithm
    class AlgorithmError : public std::runtime_error
    {
    public:
        /// @brief Exception for algorithmic errors accepting a string
        /// @param msg the error message string
        explicit AlgorithmError(const std::string& msg) : runtime_error(msg) {}

        /// @brief Exception for algorithmic errors accepting a char array
        /// @param msg the pointer to the error message char array
        explicit AlgorithmError(const char* msg) : runtime_error(msg) {}
    };

    /// @brief Custom exception to describe an error caused by an invalid mesh at a specific location
    class MeshGeometryError : public std::runtime_error
    {
    public:
        /// @brief Default exception for an invalid mesh
        MeshGeometryError() : runtime_error(""), m_invalidIndex(constants::missing::sizetValue), m_location(Mesh::Location::Unknown){};

        /// @brief Exception for an invalid mesh accepting a message string
        /// @param msg the error message string
        /// @param invalidIndex the index of the invalid mesh entity
        /// @param location the location of the invalid mesh entity (Faces, Nodes, or Edges)
        MeshGeometryError(const std::string& msg, size_t invalidIndex, Mesh::Location location) : runtime_error(msg), m_invalidIndex(invalidIndex), m_location(location)
        {
        }

        /// @brief Exception for an invalid mesh accepting a message char array pointer
        /// @param msg the pointer to the error message array
        /// @param invalidIndex the index of the invalid mesh entity
        /// @param location the location of the invalid mesh entity (Faces, Nodes, or Edges)
        MeshGeometryError(const char* msg, size_t invalidIndex, Mesh::Location location) : runtime_error(msg), m_invalidIndex(invalidIndex), m_location(location) {}

        size_t m_invalidIndex;     ///< The invalid mesh location index
        Mesh::Location m_location; ///< The location type
    };

} // namespace meshkernel
