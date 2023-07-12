//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2021.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

#pragma once

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Mesh.hpp"

#include <algorithm>
#include <exception>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

#if USE_LIBFMT
#include <fmt/format.h>
#else
#include <format>
#endif

#define STRINGIFY(str) #str
#define TO_STR_LITERAL(str) STRINGIFY(str)

namespace meshkernel
{

    /// @brief A class for generating a variadic error message.
    class VariadicErrorMessage final
    {
    public:
        /// @brief Class constructor
        /// @param[in] message The format string of the error message.
        /// @param[in] args    Arguments to be formatted.
        template <typename... Args>
        VariadicErrorMessage(std::string_view message,
                             Args const&... args)
        {
            if (sizeof...(args) == 0)
            {
                m_fmt_message = message;
            }
            else
            {
#if USE_LIBFMT
                m_fmt_message = fmt::vformat(message, fmt::make_format_args(args...));
#else
                m_fmt_message = std::vformat(message, std::make_format_args(args...));
#endif
            }
        }

        /// @brief Gets the formatted error message.
        /// @return The formatted error message.
        std::string const& GetFormatted() const { return m_fmt_message; }

    private:
        std::string m_fmt_message; ///< The formatted error message.
    };

    /// @brief A class for throwing general MeshKernel exceptions.
    class MeshKernelError : public std::exception
    {
    public:
        /// @brief Class default constructor. It is deleted to prevent constructing errors without messages.
        MeshKernelError() = delete;

        /// @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] message         The variadic error message.
        /// @param[in] source_location The source location.
        MeshKernelError(VariadicErrorMessage const& message,
                        std::source_location const& source_location = std::source_location::current())
            : m_fmt_message(message.GetFormatted()),
              m_source_location(source_location)
        {
        }

        /// @brief Class constructor parametrized by a string error message and optionally the source location.
        /// @param[in] message         The string error message.
        /// @param[in] source_location The source location.
        MeshKernelError(std::string_view message,
                        std::source_location const& source_location = std::source_location::current())
            : m_fmt_message(message),
              m_source_location(source_location)

        {
        }

        /// @brief Class destructor.
        virtual ~MeshKernelError() = default;

        /// @brief Returns the explanatory string of the error.
        /// @return the explanatory string of the error.
        const char* what() const noexcept override
        {
            std::ostringstream oss;
            oss << "Exception of type '"
                << Category()
                << "' in "
                << StrippedFilePath()
                << " ("
                << m_source_location.line()
                << ':'
                << m_source_location.column()
                << ") "
                << m_source_location.function_name()
                << ": "
                << m_fmt_message
                << '\n';
            m_fmt_message = oss.str();
            return m_fmt_message.c_str();
        }

    protected:
        /// @brief Returns the error category.
        /// @return The error category.
        virtual std::string Category() const { return "MeshKernelError"; }

    private:
        mutable std::string m_fmt_message;      ///< The formatted message
        std::source_location m_source_location; ///< The source location

        /// @brief Strips CMAKE_SRC_DIR from the path of the file name obtained from the source location.
        /// @return The stripped file name.
        std::string StrippedFilePath() const
        {
            std::string path = m_source_location.file_name();
#ifdef CMAKE_SRC_DIR
            std::string path_to_erase(TO_STR_LITERAL(CMAKE_SRC_DIR));
#ifdef _WIN32
            std::replace(path_to_erase.begin(), path_to_erase.end(), '/', '\\');
#endif
            if (size_t pos = path.find(path_to_erase); pos != std::string::npos)
            {
                // erase including the trailing slash
                path.erase(pos, path_to_erase.length() + 1);
#ifdef _WIN32
                std::replace(path.begin(), path.end(), '\\', '/');
#endif
            }
#endif
            return path;
        }
    };

    /// @brief A class for throwing "not implemented" exceptions.
    class NotImplemented final : public MeshKernelError
    {
    public:
        // @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] message         The variadic error message.
        /// @param[in] source_location The source location.
        NotImplemented(VariadicErrorMessage const& message,
                       std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location)
        {
        }

        /// @brief Class constructor parametrized by a string error message and optionally the source location.
        /// @param[in] message         The string error message.
        /// @param[in] source_location The source location.
        NotImplemented(std::string_view message,
                       std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location)
        {
        }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override { return "NotImplemented"; }
    };

    /// @brief A class for throwing algorithm exceptions
    class AlgorithmError final : public MeshKernelError
    {
    public:
        // @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] message         The variadic error message.
        /// @param[in] source_location The source location.
        AlgorithmError(VariadicErrorMessage const& message,
                       std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location)
        {
        }

        /// @brief Class constructor parametrized by a string error message and optionally the source location.
        /// @param[in] message         The string error message.
        /// @param[in] source_location The source location.
        AlgorithmError(std::string_view message,
                       std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location)
        {
        }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override { return "AlgorithmError"; }
    };

    /// @brief A class for throwing mesh geometry errors.
    class MeshGeometryError final : public MeshKernelError
    {
    public:
        // @brief Class constructor parametrized optionally by the source location.
        /// @param[in] source_location The source location.
        explicit MeshGeometryError(std::source_location const& source_location = std::source_location::current())
            : MeshKernelError("", source_location),
              m_invalid_index{constants::missing::uintValue},
              m_mesh_location(Mesh::Location::Unknown)
        {
        }

        // @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] message         The variadic error message.
        /// @param[in] invalid_index   The invalid mesh location index.
        /// @param[in] mesh_location   The location type.
        /// @param[in] source_location The source location.
        MeshGeometryError(VariadicErrorMessage const& message,
                          UInt invalid_index,
                          Mesh::Location mesh_location,
                          std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location),
              m_invalid_index(invalid_index),
              m_mesh_location(mesh_location)
        {
        }

        // @brief Class constructor parametrized by a string error message and optionally the source location.
        /// @param[in] message         The string error message.
        /// @param[in] invalid_index   The invalid mesh location index.
        /// @param[in] mesh_location   The location type.
        /// @param[in] source_location The source location.
        MeshGeometryError(std::string_view message,
                          UInt invalid_index,
                          Mesh::Location mesh_location,
                          std::source_location const& source_location = std::source_location::current())
            : MeshKernelError(message, source_location),
              m_invalid_index(invalid_index),
              m_mesh_location(mesh_location)
        {
        }

        /// @brief Returns the invalid index.
        /// @return The invalid index.
        UInt InavlidIndex() const { return m_invalid_index; }

        /// @brief Returns the mesh location.
        /// @return The mesh location.
        Mesh::Location MeshLocation() const { return m_mesh_location; }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override { return "MeshGeometryError"; }

        UInt m_invalid_index;           ///< The invalid mesh location index.
        Mesh::Location m_mesh_location; ///< The location type.
    };

} // namespace meshkernel
