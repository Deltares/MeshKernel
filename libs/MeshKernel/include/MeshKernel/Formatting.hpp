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

namespace meshkernel
{
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
            // wrirte opening braces (must be escaped)
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

} // namespace meshkernel