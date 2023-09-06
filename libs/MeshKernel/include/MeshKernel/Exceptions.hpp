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

#if USE_LIBFMT
#include <fmt/format.h>
namespace fmt_ = fmt;
#else
#include <format>
namespace fmt_ = std;
#endif

#define STRINGIFY(str) #str
#define TO_STR_LITERAL(str) STRINGIFY(str)

namespace meshkernel
{

    enum ExitCode
    {
        Success = 0,
        MeshKernelErrorCode = 1,
        NotImplementedCode = 2,
        AlgorithmErrorCode = 3,
        ConstraintErrorCode = 4,
        MeshGeometryErrorCode = 5,
        StdLibExceptionCode = 6,
        UnknownExceptionCode = 7
    };

    class ErrorCategory
    {
    public:
        ErrorCategory(std::string name,
                      ExitCode code)
            : m_name(name),
              m_code{code}
        {
        }

        std::string Name() const { return m_name; }

        ExitCode Code() const { return m_code; }

    private:
        std::string m_name;
        ExitCode m_code;
    };

    class Message final
    {
    public:
        /// @brief Class constructor
        /// @param[in] format          The format and source location
        /// @param[in] source_location The source location
        Message(const char* format,
                std::source_location const& source_location = std::source_location::current())
            : m_format(format),
              m_source_location(source_location)
        {
        }

        /// @brief Returns the format
        /// @return The format
        const char* Format() const { return m_format; }

        /// @brief Returns the source location
        /// @return The source location
        std::source_location const& SourceLocation() const { return m_source_location; }

    private:
        const char* m_format; ///< The format

        std::source_location m_source_location; ///< The source location
    };

    class MeshKernelError : public std::exception
    {
    public:
        /// @brief Class constructor
        /// @param[in] message The format and source location
        /// @param[in] args    The arguments to be formatted.
        template <typename... Args>
        MeshKernelError(Message const& message, Args&&... args)
            : m_fmt_message(),
              m_what(),
              m_source_location(message.SourceLocation())
        {
            if (sizeof...(args) == 0)
            {
                m_fmt_message = message.Format();
            }
            else
            {
                m_fmt_message = fmt_::vformat(message.Format(), fmt_::make_format_args(args...));
            }
        }

        /// @brief Class destructor.
        virtual ~MeshKernelError() = default;

        /// @brief Returns the explanatory string of the error.
        /// @return the explanatory string of the error.
        const char* what() const noexcept override
        {
            m_what = fmt_::format("Exception of type '{}' in {} ({}:{}) {}: {}\n",
                                  Category().Name(),
                                  StrippedFilePath(),
                                  m_source_location.line(),
                                  m_source_location.column(),
                                  m_source_location.function_name(),
                                  FormattedMessage());
            return m_what.c_str();
        }

        ExitCode Code() const { return Category().Code(); }

    protected:
        std::string m_fmt_message; ///< The formatted message

        /// @brief Returns the error category.
        /// @return The error category.
        virtual ErrorCategory Category() const { return {"MeshKernelError", ExitCode::MeshKernelErrorCode}; }

        /// @brief Returns the message.
        virtual std::string FormattedMessage() const { return m_fmt_message; }

    private:
        mutable std::string m_what; ///< C-style formatted message

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

    /// @brief A class for throwing algorithm exceptions
    class NotImplemented final : public MeshKernelError
    {
    public:
        // inherit
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        ErrorCategory Category() const override { return {"NotImplemented", ExitCode::NotImplementedCode}; }
    };

    /// @brief A class for throwing algorithm exceptions
    class AlgorithmError final : public MeshKernelError
    {
    public:
        // inherit
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        ErrorCategory Category() const override { return {"AlgorithmError", ExitCode::AlgorithmErrorCode}; }
    };

    /// @brief An exception class thrown when an attempt is made that violates a range constraint.
    ///
    /// 1. When an index is out of bounds or violates a range constraint.
    /// 2. Attempt to retrieve a component that does not exist.
    class ConstraintError final : public MeshKernelError
    {
    public:
        /// @brief ConstraintError constructor
        using MeshKernelError::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        ErrorCategory Category() const override { return {"ConstraintError", ExitCode::ConstraintErrorCode}; }
    };

    /// @brief A class for throwing mesh geometry errors.
    class MeshGeometryError final : public MeshKernelError
    {
    public:
        /// @brief Class constructor
        /// @param[in] invalid_index               The invalid mesh location index.
        /// @param[in] mesh_location               The location type.
        /// @param[in] message The format and source location
        /// @param[in] args                        The arguments to be formatted.
        template <typename... Args>
        MeshGeometryError(meshkernel::UInt invalid_mesh_index,
                          Mesh::Location invalid_mesh_location,
                          Message const& message,
                          Args&&... args)
            : MeshKernelError(message, std::forward<Args>(args)...),
              m_invalid_mesh_index(invalid_mesh_index),
              m_invalid_mesh_location(invalid_mesh_location)
        {
        }

        /// @brief Returns the invalid index.
        /// @return The invalid index.
        meshkernel::UInt InavlidMeshIndex() const { return m_invalid_mesh_index; }

        /// @brief Returns the mesh location.
        /// @return The mesh location.
        Mesh::Location InvalidMeshLocation() const { return m_invalid_mesh_location; }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        ErrorCategory Category() const override { return {"MeshGeometryError", ExitCode::MeshGeometryErrorCode}; }

        /// @brief Returns the message.
        std::string FormattedMessage() const override
        {
            return fmt_::format("Error occurred at index {} (location: {}). {}",
                                m_invalid_mesh_index,
                                Mesh::LocationToString.at(m_invalid_mesh_location),
                                MeshKernelError::m_fmt_message);
        }

        meshkernel::UInt m_invalid_mesh_index;  ///< The invalid mesh location index.
        Mesh::Location m_invalid_mesh_location; ///< The location type.
    };

} // namespace meshkernel