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
#include <map>
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

    /// @brief A class for throwing general MeshKernel exceptions.
    template <typename... Args>
    class MeshKernelError : public std::exception
    {
    public:
        /// @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] message         The variadic error message.
        /// @param[in] args            The arguments to be formatted.
        /// @param[in] source_location The source location.
        MeshKernelError(std::string_view message,
                        Args&&... args,
                        std::source_location const& source_location = std::source_location::current())
            : m_fmt_message(),
              m_what(),
              m_source_location(source_location)
        {
            if (sizeof...(Args) == 0)
            {
                m_fmt_message = message;
            }
            else
            {
                m_fmt_message = fmt_::vformat(message, fmt_::make_format_args(std::forward<Args>(args)...));
            }
        }

        /// @brief Class destructor.
        virtual ~MeshKernelError() = default;

        /// @brief Returns the explanatory string of the error.
        /// @return the explanatory string of the error.
        const char* what() const noexcept override
        {
            m_what = fmt_::format("Exception of type '{}' in {} ({}:{}) {}: {}\n",
                                  Category(),
                                  StrippedFilePath(),
                                  m_source_location.line(),
                                  m_source_location.column(),
                                  m_source_location.function_name(),
                                  Message());
            return m_what.c_str();
        }

    protected:
        std::string m_fmt_message; ///< The formatted message

        /// @brief Returns the error category.
        /// @return The error category.
        virtual std::string Category() const { return "MeshKernelError"; }

        /// @brief Returns the message.
        virtual std::string Message() const { return m_fmt_message; }

    private:
        mutable std::string m_what;             ///< C-style formatted message
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

    /// @brief Template deduction guide for MeshKernelError
    template <typename... Args>
    MeshKernelError(std::string_view message, Args&&... args) -> MeshKernelError<Args...>;

    /// @brief A class for throwing algorithm exceptions
    template <typename... Args>
    class NotImplemented final : public MeshKernelError<Args...>
    {
    public:
        // inherit
        using MeshKernelError<Args...>::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override
        {
            return "NotImplemented";
        }
    };

    /// @brief Template deduction guide for NotImplemented
    template <typename... Args>
    NotImplemented(std::string_view message, Args&&... args) -> NotImplemented<Args...>;

    /// @brief A class for throwing algorithm exceptions
    template <typename... Args>
    class AlgorithmError final : public MeshKernelError<Args...>
    {
    public:
        // inherit
        using MeshKernelError<Args...>::MeshKernelError;

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override
        {
            return "AlgorithmError";
        }
    };

    /// @brief Template deduction guide for AlgorithmError
    template <typename... Args>
    AlgorithmError(std::string_view message, Args&&... args) -> AlgorithmError<Args...>;

    /// @brief A class for throwing mesh geometry errors.
    template <typename... Args>
    class MeshGeometryError final : public MeshKernelError<Args...>
    {
    public:
        /// @brief Class constructor parametrized by a variadic error message and optionally the source location.
        /// @param[in] invalid_index   The invalid mesh location index.
        /// @param[in] mesh_location   The location type.
        /// @param[in] message         The variadic error message.
        /// @param[in] args            The arguments to be formatted.
        /// @param[in] source_location The source location.
        MeshGeometryError(size_t invalid_index,
                          Mesh::Location mesh_location,
                          std::string_view message,
                          Args&&... args,
                          std::source_location const& source_location = std::source_location::current())
            : MeshKernelError<Args...>(message, std::forward<Args>(args)..., source_location),
              m_invalid_index(invalid_index),
              m_mesh_location(mesh_location)
        {
        }

        /// @brief Returns the invalid index.
        /// @return The invalid index.
        size_t InavlidIndex() const { return m_invalid_index; }

        /// @brief Returns the mesh location.
        /// @return The mesh location.
        Mesh::Location MeshLocation() const { return m_mesh_location; }

    private:
        /// @brief Returns the error category.
        /// @return The  error category.
        std::string Category() const override { return "MeshGeometryError"; }

        /// @brief Returns the message.
        std::string Message() const override
        {
            std::string const message = fmt_::format("Error occurred at index {} (location: {}). {}",
                                                     m_invalid_index,
                                                     Mesh::LocationToString.at(m_mesh_location),
                                                     MeshKernelError<Args...>::m_fmt_message);
            return message;
        }

        size_t m_invalid_index;         ///< The invalid mesh location index.
        Mesh::Location m_mesh_location; ///< The location type.
    };

    /// @brief Template deduction guide for MeshGeometryError
    template <typename... Args>
    MeshGeometryError(size_t invalid_index,
                      Mesh::Location mesh_location,
                      std::string_view message,
                      Args&&... args) -> MeshGeometryError<Args...>;

} // namespace meshkernel