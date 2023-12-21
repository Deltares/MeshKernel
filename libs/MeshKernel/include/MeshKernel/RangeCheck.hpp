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

#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Formatting.hpp"

#include <algorithm>
#include <concepts>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace meshkernel
{
    namespace range_check
    {
        /// @brief Defines the checkable types: floating point types and inetgral types except for bool and char.
        template <typename T>
        concept RangeCheckableType = std::floating_point<T> ||
                                     (std::integral<T> &&
                                      !std::same_as<T, bool> &&
                                      !std::same_as<T, char> &&
                                      !std::same_as<T, char8_t> &&
                                      !std::same_as<T, char16_t> &&
                                      !std::same_as<T, char32_t> &&
                                      !std::same_as<T, wchar_t>);

        /// @brief  Keys of performed comparison
        enum class Comparison
        {
            Equal,                   ///< Equal to
            NotEqual,                ///< Not equal to
            Greater,                 ///< Greater than
            GreaterEqual,            ///< Greater than or equal to
            Less,                    ///< Less than
            LessEqual,               ///< Less than or equal to
            InClosedInterval,        ///< In closed interval
            InOpenInterval,          ///< In open open interval
            InRightHalfOpenInterval, ///< In semi-open from above interval
            InLeftHalfOpenInterval,  ///< In semi-open from below interval
            OutsideClosedInterval,   ///< Outside closed interval
            OutsideOpenInterval,     ///< Outside open interval
            OneOf,                   ///< One of
            NoneOf                   ///< None of
        };

        /// @brief Maps the comparison keys to valid range format string (used for the generation of error messages)
        inline static std::unordered_map<Comparison, std::string> const ValidRangeFormat = {
            {Comparison::Equal, "value = {}"},
            {Comparison::NotEqual, "value != {}"},
            {Comparison::Greater, "value > {}"},
            {Comparison::GreaterEqual, "value >= {}"},
            {Comparison::Less, "value < {}"},
            {Comparison::LessEqual, "value <= {}"},
            {Comparison::InClosedInterval, "{} <= value <= {}"},
            {Comparison::InOpenInterval, "{} < value < {}"},
            {Comparison::InRightHalfOpenInterval, "{} <= value < {}"},
            {Comparison::InLeftHalfOpenInterval, "{} < value <= {}"},
            {Comparison::OutsideClosedInterval, "value < {} and value > {}"},
            {Comparison::OutsideOpenInterval, "value <= {} and value >= {}"},
            {Comparison::OneOf, "value is one of {}"},
            {Comparison::NoneOf, "value is none of {}"} //
        };

        /// @brief Checks the validity of a value given a bound, supports the predicates ==, !=, >, >=, <, and <=
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param bound The bound against which the value is compared
        /// @param predicate The comparison predicate (one of ==, !=, >, >=, <, and <=)
        /// @param comparison The key of the performed comparison
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckRange(T const& value,
                                      T const& bound,
                                      std::function<bool(T const&, T const&)> predicate,
                                      Comparison const comparison,
                                      std::string_view const variable_name)
        {
            if (!predicate(value, bound))
            {
                throw RangeError(
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.",
                                   ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    bound);
            }
        }

        /// @brief Checks the validity of a value given interval bounds, supports the predicates ==, !=, >, >=, <, and <=
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param predicate The comparison predicate (one of <= ... <=, < ... <, < ... <=, and <= ... <)
        /// @param comparison The key of the performed comparison
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckRange(T const& value,
                                      T const& interval_lower_bound,
                                      T const& interval_upper_bound,
                                      std::function<bool(T const&, T const&, T const&)> predicate,
                                      Comparison const comparison,
                                      std::string_view const variable_name)
        {
            if (!predicate(value, interval_lower_bound, interval_upper_bound))
            {
                throw RangeError(
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.",
                                   ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    interval_lower_bound,
                    interval_upper_bound);
            }
        }

        /// @brief Checks the validity of a value given a vector of values, supports the predicates one of and none of
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param values The vector of values
        /// @param predicate The comparison predicate (one of <= ... <=, < ... <, < ... <=, and <= ... <)
        /// @param comparison The key of the performed comparison
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckRange(T const& value,
                                      std::vector<T> const& values,
                                      std::function<bool(T const&, std::vector<T> const&)> predicate,
                                      Comparison const comparison,
                                      std::string_view const variable_name)
        {
            if (!predicate(value, values))
            {
                throw RangeError(
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.",
                                   ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    values);
            }
        }

        /// @brief Checks if two values are equal
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param other_value The other value
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckEqual(T const& value,
                                      T const& other_value,
                                      std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          other_value,
                          std::equal_to<T>(),
                          Comparison::Equal,
                          variable_name);
        }

        /// @brief Checks if two values are not equal
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param other_value The other value
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckNotEqual(T const& value,
                                         T const& other_value,
                                         std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          other_value,
                          std::not_equal_to<T>(),
                          Comparison::NotEqual,
                          variable_name);
        }

        /// @brief Checks if a value is strictly greater than another
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param bound The bound
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckGreater(T const& value,
                                        T const& bound,
                                        std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::greater<T>(),
                          Comparison::Greater,
                          variable_name);
        }

        /// @brief Checks if a value is greater than or equal to another
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param bound The bound
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckGreaterEqual(T const& value,
                                             T const& bound,
                                             std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::greater_equal<T>(),
                          Comparison::GreaterEqual,
                          variable_name);
        }

        /// @brief Checks if a value is strictly less than another
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param bound The bound
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckLess(T const& value,
                                     T const& bound,
                                     std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::less<T>(),
                          Comparison::Less,
                          variable_name);
        }

        /// @brief Checks if a value is less than or equal to another
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param bound The bound
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckLessEqual(T const& value,
                                          T const& bound,
                                          std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::less_equal<T>(),
                          Comparison::LessEqual,
                          variable_name);
        }

        /// @brief Checks if a value is in a closed inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInClosedInterval(T const& value,
                                                 T const& interval_lower_bound,
                                                 T const& interval_upper_bound,
                                                 std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less_equal{}(l, v) && std::less_equal{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::InClosedInterval,
                          variable_name);
        }

        /// @brief Checks if a value is in a closed inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param inteval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInClosedInterval(T const& value,
                                                 std::pair<T, T> const& inteval_bounds,
                                                 std::string_view const variable_name)
        {
            CheckInClosedInterval(value,
                                  inteval_bounds.first,
                                  inteval_bounds.second,
                                  variable_name);
        }

        /// @brief Checks if a value is in an open inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInOpenInterval(T const& value,
                                               T const& interval_lower_bound,
                                               T const& interval_upper_bound,
                                               std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less{}(l, v) && std::less{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::InOpenInterval,
                          variable_name);
        }

        /// @brief Checks if a value is in an open inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param inteval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInOpenInterval(T const& value,
                                               std::pair<T, T> const& interval_bounds,
                                               std::string_view const variable_name)
        {
            CheckInOpenInterval(value,
                                interval_bounds.first,
                                interval_bounds.second,
                                variable_name);
        }

        /// @brief Checks if a value is in an open from above and closed from below inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInRightHalfOpenInterval(T const& value,
                                                        T const& interval_lower_bound,
                                                        T const& interval_upper_bound,
                                                        std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less_equal{}(l, v) && std::less{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::InRightHalfOpenInterval,
                          variable_name);
        }

        /// @brief Checks if a value is in an open from above and closed from below inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInRightHalfOpenInterval(T const& value,
                                                        std::pair<T, T> const& interval_bounds,
                                                        std::string_view const variable_name)
        {
            CheckInRightHalfOpenInterval(value,
                                         interval_bounds.first,
                                         interval_bounds.second,
                                         variable_name);
        }

        /// @brief Checks if a value is in an open from below and closed from above inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInLeftHalfOpenInterval(T const& value,
                                                       T const& interval_lower_bound,
                                                       T const& interval_upper_bound,
                                                       std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less{}(l, v) && std::less_equal{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::InLeftHalfOpenInterval,
                          variable_name);
        }

        /// @brief Checks if a value is in an open from below and closed from above inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckInLeftHalfOpenInterval(T const& value,
                                                       std::pair<T, T> const& interval_bounds,
                                                       std::string_view const variable_name)
        {
            CheckInLeftHalfOpenInterval(value,
                                        interval_bounds.first,
                                        interval_bounds.second,
                                        variable_name);
        }

        /// @brief Checks if a value is outsidea closed inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckOutsideClosedInterval(T const& value,
                                                      T const& interval_lower_bound,
                                                      T const& interval_upper_bound,
                                                      std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less{}(v, l) || std::greater{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::OutsideClosedInterval,
                          variable_name);
        }

        /// @brief Checks if a value is outsidea closed inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckOutsideClosedInterval(T const& value,
                                                      std::pair<T, T> const& interval_bounds,
                                                      std::string_view const variable_name)
        {
            CheckOutsideClosedInterval(value,
                                       interval_bounds.first,
                                       interval_bounds.second,
                                       variable_name);
        }

        /// @brief Checks if a value is outsidean open inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_lower_bound The lower bound of the interval
        /// @param interval_upper_bound The upper bound of the interval
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckOutsideOpenInterval(T const& value,
                                                    T const& interval_lower_bound,
                                                    T const& interval_upper_bound,
                                                    std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less_equal{}(v, l) || std::greater_equal{}(v, u);
            };

            CheckRange<T>(value,
                          interval_lower_bound,
                          interval_upper_bound,
                          predicate,
                          Comparison::OutsideClosedInterval,
                          variable_name);
        }

        /// @brief Checks if a value is outsidean open inetrval
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param interval_bounds The interval bounds
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckOutsideOpenInterval(T const& value,
                                                    std::pair<T, T> const& interval_bounds,
                                                    std::string_view const variable_name)
        {
            CheckOutsideOpenInterval(value,
                                     interval_bounds.first,
                                     interval_bounds.second,
                                     variable_name);
        }

        /// @brief Checks if a value is an element of a vector of values
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param values The vector of values
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckOneOf(T const& value,
                                      std::vector<T> const& values,
                                      std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, std::vector<T> const& l)
            {
                return std::any_of(l.begin(), l.end(), [&v](T const& i)
                                   { return std::equal_to<T>{}(i, v); });
            };

            CheckRange<T>(value,
                          values,
                          predicate,
                          Comparison::OneOf,
                          variable_name);
        }

        /// @brief Checks if a value is not an element of a vector of values
        /// @tparam T The type of the value to check (RangeCheckableType)
        /// @param value The value to check
        /// @param values The vector of values
        /// @param variable_name The name or description of the variable whose value is checked
        ///                      (It is recommended to use a concise string)
        template <RangeCheckableType T>
        inline static void CheckNoneOf(T const& value,
                                       std::vector<T> const& values,
                                       std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, std::vector<T> const& l)
            {
                return std::none_of(l.begin(), l.end(), [&v](T const& i)
                                    { return std::equal_to<T>{}(i, v); });
            };

            CheckRange<T>(value,
                          values,
                          predicate,
                          Comparison::NoneOf,
                          variable_name);
        }

    } // namespace range_check

} // namespace meshkernel
