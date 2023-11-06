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

#include <cmath>

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Utilities/NumericFunctions.hpp"
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{

    /// @brief A struct describing a point in a two-dimensional space
    class Point
    {
    public:
        double x; ///< X-coordinate
        double y; ///< Y-coordinate

        /// @brief Constructor initializing with missing values
        Point()
            : x(constants::missing::doubleValue),
              y(constants::missing::doubleValue)
        {
        }

        /// @brief Constructor initializing with given coordinates
        /// @param[in] x coordinate
        /// @param[in] y coordinate
        Point(double x, double y)
            : x(x),
              y(y)
        {
        }

        /// @brief Inplace add point to point.
        Point& operator+=(const Point& p);

        /// @brief Inplace subtract point from point.
        Point& operator-=(const Point& p);

        /// @brief Inplace add vector to a point.
        Point& operator+=(const Vector& vec);

        /// @brief Inplace subtract vector from a point.
        Point& operator-=(const Vector& vec);

        /// @brief Inplace divide point by point.
        ///
        /// \note No check is made that any component of the divisor is 0
        Point& operator/=(const Point& p);

        /// @brief Inplace multiply point by point.
        Point& operator*=(const Point& p);

        /// @brief Inplace add scalar to point.
        Point& operator+=(const double p);

        /// @brief Inplace subtract scalar from point.
        Point& operator-=(const double p);

        /// @brief Inplace divide point by scalar.
        ///
        /// \note No check is made that the divisor is 0
        Point& operator/=(const double p);

        /// @brief Inplace multiply point by scalar.
        Point& operator*=(const double p);

        /// @brief Overloads multiplication with a integer
        [[nodiscard]] Point operator*(int const& rhs) const { return Point(x * rhs, y * rhs); }

        /// @brief Transforms spherical coordinates to cartesian
        void TransformSphericalToCartesian(double referenceLatitude)
        {
            static double constexpr m_trans_factor =
                constants::conversion::degToRad *
                constants::geometric::earth_radius; ///< Factor used in the transformation from spherical to Cartesian coordinates

            x = x * m_trans_factor * std::cos(constants::conversion::degToRad * referenceLatitude);
            y = y * m_trans_factor;
        }

        /// @brief Determines if one of the point coordinates equals to \p missingValue
        [[nodiscard]] bool IsValid(const double missingValue = constants::missing::doubleValue) const
        {
            const bool isInvalid = IsEqual(x, missingValue) ||
                                   IsEqual(y, missingValue) ||
                                   IsEqual(x, constants::missing::innerOuterSeparator) ||
                                   IsEqual(y, constants::missing::innerOuterSeparator);

            return !isInvalid;
        }

        /// @brief Set the point to be invalid.
        void SetInvalid();
    };

    /// @brief Unary minus
    ///
    /// @returns \f$ (-p1.x, -p1.y)\f$
    Point operator-(const Point& pnt);

    /// @brief Add two points
    ///
    /// @returns \f$ (p1.x + p2.x, p1.y + p2.y)\f$
    Point operator+(const Point& p1, const Point& p2);

    /// @brief Add vector to a point
    ///
    /// @returns \f$ (p.x + v.x, p.y + v.y)\f$
    Point operator+(const Point& pnt, const Vector& vec);

    /// @brief Add points and scalar
    ///
    /// @returns \f$ (p.x + x, p.y + x)\f$
    Point operator+(const Point& p, const double x);

    /// @brief Add points and scalar
    ///
    /// @returns \f$ (p.x + x, p.y + x)\f$
    Point operator+(const double x, const Point& p);

    /// @brief Subtract two points
    ///
    /// @returns \f$ (p1.x - p2.x, p1.y - p2.y)\f$
    Point operator-(const Point& p1, const Point& p2);

    /// @brief Subtract vector from a point
    ///
    /// @returns \f$ (p.x - v.x, p.y - v.y)\f$
    Point operator-(const Point& pnt, const Vector& vec);

    /// @brief Subtract scalar from point
    ///
    /// @returns \f$ (p.x - x, p.y - x)\f$
    Point operator-(const Point& p, const double x);

    /// @brief Subtract point from scalar
    ///
    /// @returns \f$ (x - p.x, x - p.y)\f$
    Point operator-(const double x, const Point& p);

    /// @brief Multiply point by a point
    ///
    /// Computes a point-wise product.
    /// @returns \f$ (p1.x * p2.x, p1.y * p2.y)\f$
    Point operator*(const Point& p1, const Point& p2);

    /// @brief Multiply point by a scalar
    /// @returns \f$ (p.x * x, p.y * x)\f$
    Point operator*(const double x, const Point& p);

    /// @brief Multiply point by a scalar
    /// @returns \f$ (p.x * x, p.y * x)\f$
    Point operator*(const Point& p, const double x);

    /// @brief Divide point by a point
    ///
    /// Computes a point-wise division.
    /// @returns \f$ (p1.x / p2.x, p1.y / p2.y)\f$
    /// \note No check is made to determine is divisors are zero.
    Point operator/(const Point& p1, const Point& p2);

    /// @brief Divide point by a scalar
    /// \note No check is made to determine is divisor is zero.
    /// @returns \f$ (p.x / x, p.y / x)\f$
    Point operator/(const Point& p, const double x);

    /// @brief Test points for equality using a default tolerance
    /// @returns \f$ p1.x = p2.x \wedge p1.y = p2.y)\f$
    bool operator==(const Point& p1, const Point& p2);

    /// @brief Test points for equality using a default tolerance
    /// @returns \f$ p1.x = p2.x \wedge p1.y = p2.y)\f$
    bool operator!=(const Point& p1, const Point& p2);

    /// @brief Test points for equality upto a tolerance
    /// @param [in] p1 First point to compare
    /// @param [in] p2 Second point to compare
    /// @param[in] epsilon Relative tolerance to which the values are compared.
    /// @returns Boolean value indicating where the points p1 and p2 are equal upto a relative tolerance.
    bool IsEqual(const Point& p1, const Point& p2, const double epsilon);

    /// @brief Compute the point at some position along the line connecting start- and end-point.
    ///
    /// \f$ r = (1-\lambda) s + \lambda e \f$
    Point PointAlongLine(const Point& startPoint, const Point& endPoint, const double lambda);

    /// @brief Rotate a point around a reference
    /// @param[in] point The point to rotate
    /// @param[in] angle The rotation angle
    /// @param[in] reference The reference point where rotation should be performed
    /// @returns The rotated point
    static Point Rotate(const Point& point, const double angle, const Point& reference)
    {
        const auto translatedPoint = point - reference;

        const auto angleInRad = angle * constants::conversion::degToRad;
        const auto cosineAngle = std::cos(angleInRad);
        const auto sinAngle = std::sin(angleInRad);
        Point result(translatedPoint.x * cosineAngle - translatedPoint.y * sinAngle,
                     translatedPoint.x * sinAngle + translatedPoint.y * cosineAngle);

        return result + reference;
    }

    /// @brief Get the delta-x in Cartesian coordinate system
    double GetDeltaXCartesian(const Point& p1, const Point& p2);

    /// @brief Get the delta-y in Cartesian coordinate system
    double GetDeltaYCartesian(const Point& p1, const Point& p2);

    /// @brief Get the delta-x and -y in Cartesian coordinate system
    Vector GetDeltaCartesian(const Point& p1, const Point& p2);

} // namespace meshkernel

inline meshkernel::Point& meshkernel::Point::operator+=(const Point& p)
{
    x += p.x;
    y += p.y;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator-=(const Point& p)
{
    x -= p.x;
    y -= p.y;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator+=(const Vector& vec)
{
    x += vec.x();
    y += vec.y();
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator-=(const Vector& vec)
{
    x -= vec.x();
    y -= vec.y();
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator/=(const Point& p)
{
    x /= p.x;
    y /= p.y;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator*=(const Point& p)
{
    x *= p.x;
    y *= p.y;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator+=(const double p)
{
    x += p;
    y += p;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator-=(const double p)
{
    x -= p;
    y -= p;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator/=(const double p)
{
    x /= p;
    y /= p;
    return *this;
}

inline meshkernel::Point& meshkernel::Point::operator*=(const double p)
{
    x *= p;
    y *= p;
    return *this;
}

inline meshkernel::Point meshkernel::operator-(const Point& pnt)
{
    return Point(-pnt.x, -pnt.y);
}

inline meshkernel::Point meshkernel::operator+(const Point& p1, const Point& p2)
{
    return Point(p1.x + p2.x, p1.y + p2.y);
}

inline meshkernel::Point meshkernel::operator+(const Point& pnt, const Vector& vec)
{
    return Point(pnt.x + vec.x(), pnt.y + vec.y());
}

inline meshkernel::Point meshkernel::operator+(const Point& p, const double x)
{
    return Point(p.x + x, p.y + x);
}

inline meshkernel::Point meshkernel::operator+(const double x, const Point& p)
{
    return Point(p.x + x, p.y + x);
}

inline meshkernel::Point meshkernel::operator-(const Point& p1, const Point& p2)
{
    return Point(p1.x - p2.x, p1.y - p2.y);
}

inline meshkernel::Point meshkernel::operator-(const Point& pnt, const Vector& vec)
{
    return Point(pnt.x - vec.x(), pnt.y - vec.y());
}

inline meshkernel::Point meshkernel::operator-(const Point& p, const double x)
{
    return Point(p.x - x, p.y - x);
}

inline meshkernel::Point meshkernel::operator-(const double x, const Point& p)
{
    return Point(x - p.x, x - p.y);
}

inline meshkernel::Point meshkernel::operator*(const Point& p1, const Point& p2)
{
    return Point(p1.x * p2.x, p1.y * p2.y);
}

inline meshkernel::Point meshkernel::operator*(const double x, const Point& p)
{
    return Point(x * p.x, x * p.y);
}

inline meshkernel::Point meshkernel::operator*(const Point& p, const double x)
{
    return Point(x * p.x, x * p.y);
}

inline meshkernel::Point meshkernel::operator/(const Point& p1, const Point& p2)
{
    return Point(p1.x / p2.x, p1.y / p2.y);
}

inline meshkernel::Point meshkernel::operator/(const Point& p, const double x)
{
    return Point(p.x / x, p.y / x);
}

inline bool meshkernel::operator==(const Point& p1, const Point& p2)
{
    return IsEqual(p1.x, p2.x) && IsEqual(p1.y, p2.y);
}

inline bool meshkernel::operator!=(const Point& p1, const Point& p2)
{
    return !(p1 == p2);
}

inline bool meshkernel::IsEqual(const Point& p1, const Point& p2, const double epsilon)
{
    return IsEqual(p1.x, p2.x, epsilon) && IsEqual(p1.y, p2.y, epsilon);
}

inline double meshkernel::GetDeltaXCartesian(const Point& p1, const Point& p2)
{
    return p2.x - p1.x;
}

inline double meshkernel::GetDeltaYCartesian(const Point& p1, const Point& p2)
{
    return p2.y - p1.y;
}

inline meshkernel::Vector meshkernel::GetDeltaCartesian(const Point& p1, const Point& p2)
{
    return Vector(GetDeltaXCartesian(p1, p2), GetDeltaYCartesian(p1, p2));
}

inline meshkernel::Point meshkernel::PointAlongLine(const Point& startPoint, const Point& endPoint, const double lambda)
{
    return (1.0 - lambda) * startPoint + lambda * endPoint;
}

void inline meshkernel::Point::SetInvalid()
{
    x = constants::missing::doubleValue;
    y = constants::missing::doubleValue;
}
