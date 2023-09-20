//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <concepts>
#include <string>
#include <vector>

#if USE_LIBFMT
#include <fmt/format.h>
namespace fmt_ns = fmt;
#else
#include <format>
namespace fmt_ns = std;
#endif

/// @brief Defines the formattable types
template <typename T>
concept FormattableType = std::floating_point<T> ||
                          std::integral<T> ||
                          std::same_as<T, std::string>;

/// @brief Specialization of std::formatter for std::vector
/// @tparam T Vector type (must be FormattableType)
template <FormattableType T>
struct fmt_ns::formatter<std::vector<T>> : fmt_ns::formatter<T>
{
    /// @brief
    /// @tparam FormatContext
    /// @param vector The vector to be formatted
    /// @param format_context The format context
    /// @return Iterator past the end of the context's output buffer
    template <typename FormatContext>
    auto format(std::vector<T> const& vector,
                FormatContext& format_context) const
    {
        // get ref to the output buffer
        auto&& out = format_context.out();
        // write opening braces (must be escaped)
        fmt_ns::format_to(out, "{{");
        bool first = true;
        for (auto const& item : vector)
        {
            if (first)
            {
                first = false;
            }
            else
            {
                fmt_ns::format_to(out, ", ");
            }
            fmt_ns::formatter<T>::format(item, format_context);
        }
        // write closing braces and return
        return fmt_ns::format_to(out, "}}");
    }
};
