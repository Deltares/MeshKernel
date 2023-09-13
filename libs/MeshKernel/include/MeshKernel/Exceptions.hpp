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
#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh.hpp"

#include <algorithm>
#include <exception>
#include <source_location>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>

// set the namespace of the formatting library
#if USE_LIBFMT
#include <fmt/format.h>
namespace fmt_ns = fmt;
#else
#include <format>
namespace fmt_ns = std;
#endif

#define STRINGIFY(str) #str
#define TO_STR_LITERAL(str) STRINGIFY(str)

namespace meshkernel
{

    /// @brief Enumeration of exit codes
    enum ExitCode
    {
        Success = 0,                ///< Success
        MeshKernelErrorCode = 1,    ///< MehKernel error
        NotImplementedCode = 2,     ///< Not implemented error
        AlgorithmErrorCode = 3,     ///< Algorithm error
        ConstraintErrorCode = 4,    ///< Constraint error
        MeshGeometryErrorCode = 5,  ///< Geometry error
        LinearAlgebraErrorCode = 6, ///< Lienar algebra error
        StdLibExceptionCode = 7,    ///< Standrad library exception
        UnknownExceptionCode = 8    ///< Unknown exception
    };

    /// @brief Contains error category information
    class ErrorCategory
    {
    public:
        /// @brief Class constructor
        /// @param[in] name      The name of the error category
        /// @param[in] exit_code The exit code of the error category
        ErrorCategory(std::string_view name,
                      ExitCode exit_code)
            : m_name(name),
              m_exit_code{exit_code}
        {
        }

        /// @brief Returns the name of the error category
        /// @return The name of the error category
        [[nodiscard]] std::string_view Name() const { return m_name; }

        /// @brief Return the exit code of the error category
        /// @return The exit code of the error category
        [[nodiscard]] ExitCode Code() const { return m_exit_code; }

    private:
        std::string_view m_name; ///< The name of the category
        ExitCode m_exit_code;    ///< The exit code of the category
    };

    /// @brief Manages the format string and source location
    class FormatString final
    {
    public:
        /// @brief Class constructor
        /// @param[in] format_string   The format string
        /// @param[in] source_location The source location
        FormatString(const char* const format_string,
                     std::source_location const& source_location = std::source_location::current())
            : m_format_string(format_string),
              m_source_location(source_location)
        {
        }

        /// @brief Returns the format string
        /// @return The format string
        [[nodiscard]] std::string_view String() const { return m_format_string; }

        /// @brief Returns the source location
        /// @return The source location
        [[nodiscard]] std::source_location const& SourceLocation() const { return m_source_location; }

    private:
        std::string_view m_format_string; ///< The format

        std::source_location m_source_location; ///< The source location
    };

    /// @brief A class for throwing general MeshKernel exceptions
    class MeshKernelError : public std::exception
    {
    public:
        /// @brief Class constructor
        /// @param[in] format_string The format string
        /// @param[in] args         The arguments to be formatted.
        template <typename... Args>
        MeshKernelError(FormatString const& format_string, Args&&... args)
            : m_what(),
              m_source_location(format_string.SourceLocation())
        {
            if (sizeof...(args) == 0)
            {
                m_formatted_message = format_string.String();
            }
            else
            {
                m_formatted_message = fmt_ns::vformat(format_string.String(), fmt_ns::make_format_args(args...));
            }
        }

        /// @brief Class destructor.
        virtual ~MeshKernelError() = default;

        /// @brief Returns the explanatory string of the error.
        /// @return The explanatory string of the error.
        const char* what() const noexcept override
        {
#if HAVE_SRC_LOC_IN_ERR_MSGS
            m_what = fmt_ns::format("Exception of type '{}' in {} ({}:{}) {}: {}\n",
                                    Category().Name(),
                                    StrippedFilePath(),
                                    m_source_location.line(),
                                    m_source_location.column(),
                                    m_source_location.function_name(),
                                    FormattedMessage());
#else
            m_what = fmt_ns::format("Exception of type '{}': {}\n",
                                    Category().Name(),
                                    FormattedMessage());
#endif
            return m_what.c_str();
        }

        /// @brief Returns the exit code
        /// @return The exit code
        [[nodiscard]] ExitCode Code() const { return Category().Code(); }

    protected:
        /// @brief Returns the error category.
        /// @return The error category.
        [[nodiscard]] virtual ErrorCategory Category() const
        {
            return {"MeshKernelError", ExitCode::MeshKernelErrorCode};
        }

        /// @brief Returns the message.
        [[nodiscard]] virtual std::string FormattedMessage() const { return m_formatted_message; }

    private:
        std::string m_formatted_message; ///< The formatted message

        mutable std::string m_what; ///< C-style formatted message

        std::source_location m_source_location; ///< The source location

        /// @brief Strips CMAKE_SRC_DIR from the path of the file name obtained from the source location.
        /// @return The stripped file name.
        [[nodiscard]] std::string StrippedFilePath() const
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

    /// @brief A class for throwing not implemented exceptions
    class NotImplemented final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        [[nodiscard]] ErrorCategory Category() const override
        {
            return {"NotImplemented", ExitCode::NotImplementedCode};
        }
    };

    /// @brief A class for throwing algorithm exceptions
    class AlgorithmError final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        [[nodiscard]] ErrorCategory Category() const override
        {
            return {"AlgorithmError", ExitCode::AlgorithmErrorCode};
        }
    };

    /// @brief An exception class thrown when an attempt is made that violates a range constraint.
    ///
    /// 1. When an index is out of bounds or violates a range constraint.
    /// 2. Attempt to retrieve a component that does not exist.
    class ConstraintError final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        [[nodiscard]] ErrorCategory Category() const override
        {
            return {"ConstraintError", ExitCode::ConstraintErrorCode};
        }
    };

    /// @brief A class for throwing mesh geometry errors.
    class MeshGeometryError final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        /// @param[in] mesh_index    The mesh index
        /// @param[in] mesh_location The mesh location
        /// @param[in] format_string The format string
        /// @param[in] args          The arguments to be formatted.
        template <typename... Args>
        MeshGeometryError(meshkernel::UInt mesh_index,
                          Mesh::Location mesh_location,
                          FormatString const& format_string,
                          Args&&... args)
            : MeshKernelError(format_string, std::forward<Args>(args)...),
              m_mesh_index(mesh_index),
              m_mesh_location(mesh_location)
        {
        }

        /// @brief Returns the invalid index.
        /// @return The invalid index.
        [[nodiscard]] meshkernel::UInt MeshIndex() const { return m_mesh_index; }

        /// @brief Returns the mesh location.
        /// @return The mesh location.
        [[nodiscard]] Mesh::Location MeshLocation() const { return m_mesh_location; }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        [[nodiscard]] ErrorCategory Category() const override
        {
            return {"MeshGeometryError", ExitCode::MeshGeometryErrorCode};
        }

        /// @brief Returns the message.
        [[nodiscard]] std::string FormattedMessage() const override
        {
            return fmt_ns::format("Error occurred at index {} (location: {}). {}",
                                  m_mesh_index,
                                  Mesh::LocationToString.at(m_mesh_location),
                                  MeshKernelError::FormattedMessage());
        }

        meshkernel::UInt m_mesh_index;  ///< The invalid mesh location index.
        Mesh::Location m_mesh_location; ///< The location type.
    };

    /// @brief A class for throwing linear algebra exceptions
    class LinearAlgebraError final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        [[nodiscard]] ErrorCategory Category() const override
        {
            return {"LinearAlgebraError", ExitCode::LinearAlgebraErrorCode};
        }
    };

} // namespace meshkernel
