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
        template <typename T>
        concept RangeCheckableType = std::floating_point<T> ||
                                     (std::integral<T> &&
                                      !std::same_as<T, bool> &&
                                      !std::same_as<T, char> &&
                                      !std::same_as<T, char8_t> &&
                                      !std::same_as<T, char16_t> &&
                                      !std::same_as<T, char32_t> &&
                                      !std::same_as<T, wchar_t>);

        enum class Comparison
        {
            Equal,
            NotEqual,
            Greater,
            GreaterEqual,
            Less,
            LessEqual,
            InClosedInterval,
            InOpenInterval,
            InSemiOpenFromAboveInterval,
            InSemiOpenFromBelowInterval,
            OneOf,
            NoneOf
        };

        inline static std::unordered_map<Comparison, std::string> const ValidRangeFormat = {
            {Comparison::Equal, "{} = {}"},
            {Comparison::NotEqual, "{} != {}"},
            {Comparison::Greater, "{} > {}"},
            {Comparison::GreaterEqual, "{} >= {}"},
            {Comparison::Less, "{} < {}"},
            {Comparison::LessEqual, "{} <= {}"},
            {Comparison::InClosedInterval, "{} <= {} <= {}"},
            {Comparison::InOpenInterval, "{} <= {} <= {}"},
            {Comparison::InSemiOpenFromAboveInterval, "{} <= {} < {}"},
            {Comparison::InSemiOpenFromBelowInterval, "{} < {} <= {}"},
            {Comparison::OneOf, "{} is one of {}"},
            {Comparison::NoneOf, "{} is none of {}"} //
        };

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
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.", ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    variable_name,
                    bound);
            }
        }

        template <RangeCheckableType T>
        inline static void CheckRange(T const& value,
                                      T const& lower_bound,
                                      T const& upper_bound,
                                      std::function<bool(T const&, T const&, T const&)> predicate,
                                      Comparison const comparison,
                                      std::string_view const variable_name)
        {
            if (!predicate(value, lower_bound, upper_bound))
            {
                throw RangeError(
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.", ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    lower_bound,
                    variable_name,
                    upper_bound);
            }
        }

        template <RangeCheckableType T>
        inline static void CheckRange(T const& value,
                                      std::vector<T> const& list,
                                      std::function<bool(T const&, std::vector<T> const&)> predicate,
                                      Comparison const comparison,
                                      std::string_view const variable_name)
        {
            if (!predicate(value, list))
            {
                throw RangeError(
                    fmt_ns::format("{{}} = {{}} is invalid. Valid range: {}.",
                                   ValidRangeFormat.at(comparison)),
                    variable_name,
                    value,
                    variable_name,
                    list);
            }
        }

        template <RangeCheckableType T>
        inline static void CheckEqual(T const& value,
                                      T const& bound,
                                      std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::equal_to<T>(),
                          Comparison::Equal,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckNotEqual(T const& value,
                                         T const& bound,
                                         std::string_view const variable_name)
        {
            CheckRange<T>(value,
                          bound,
                          std::not_equal_to<T>(),
                          Comparison::NotEqual,
                          variable_name);
        }

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

        template <RangeCheckableType T>
        inline static void CheckInClosedInterval(T const& value,
                                                 T const& lower_bound,
                                                 T const& upper_bound,
                                                 std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less_equal{}(l, v) && std::less_equal{}(v, u);
            };

            CheckRange<T>(value,
                          lower_bound,
                          upper_bound,
                          predicate,
                          Comparison::InClosedInterval,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInClosedInterval(T const& value,
                                                 std::pair<T, T> const& bounds,
                                                 std::string_view const variable_name)
        {
            CheckInClosedInterval(value,
                                  bounds.first,
                                  bounds.second,
                                  variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInOpenInterval(T const& value,
                                               T const& lower_bound,
                                               T const& upper_bound,
                                               std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less{}(l, v) && std::less{}(v, u);
            };

            CheckRange<T>(value,
                          lower_bound,
                          upper_bound,
                          predicate,
                          Comparison::InOpenInterval,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInOpenInterval(T const& value,
                                               std::pair<T, T> const& bounds,
                                               std::string_view const variable_name)
        {
            CheckInOpenInterval(value,
                                bounds.first,
                                bounds.second,
                                variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInSemiOpenFromAboveInterval(T const& value,
                                                            T const& lower_bound,
                                                            T const& upper_bound,
                                                            std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less_equal{}(l, v) && std::less{}(v, u);
            };

            CheckRange<T>(value,
                          lower_bound,
                          upper_bound,
                          predicate,
                          Comparison::InSemiOpenFromAboveInterval,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInSemiOpenFromAboveInterval(T const& value,
                                                            std::pair<T, T> const& bounds,
                                                            std::string_view const variable_name)
        {
            CheckInSemiOpenFromAboveInterval(value,
                                             bounds.first,
                                             bounds.second,
                                             variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInSemiOpenFromBelowInterval(T const& value,
                                                            T const& lower_bound,
                                                            T const& upper_bound,
                                                            std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, T const& l, T const& u)
            {
                return std::less{}(l, v) && std::less_equal{}(v, u);
            };

            CheckRange<T>(value,
                          lower_bound,
                          upper_bound,
                          predicate,
                          Comparison::InSemiOpenFromBelowInterval,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckInSemiOpenFromBelowInterval(T const& value,
                                                            std::pair<T, T> const& bounds,
                                                            std::string_view const variable_name)
        {
            CheckInSemiOpenFromBelowInterval(value,
                                             bounds.first,
                                             bounds.second,
                                             variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckOneOf(T const& value,
                                      std::vector<T> const& list,
                                      std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, std::vector<T> const& l)
            {
                return std::any_of(l.begin(), l.end(), [&v](T const& i)
                                   { return std::equal_to<T>{}(i, v); });
            };

            CheckRange<T>(value,
                          list,
                          predicate,
                          Comparison::OneOf,
                          variable_name);
        }

        template <RangeCheckableType T>
        inline static void CheckNoneOf(T const& value,
                                       std::vector<T> const& list,
                                       std::string_view const variable_name)
        {
            auto const predicate = [](T const& v, std::vector<T> const& l)
            {
                return std::none_of(l.begin(), l.end(), [&v](T const& i)
                                    { return std::equal_to<T>{}(i, v); });
            };

            CheckRange<T>(value,
                          list,
                          predicate,
                          Comparison::NoneOf,
                          variable_name);
        }

    } // namespace range_check

} // namespace meshkernel