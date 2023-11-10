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

#include <cmath>

#include "MeshKernel/Constants.hpp"

namespace meshkernel
{
    /// @brief A class defining a vector
    class Vector
    {
    public:
        /// @brief Default constructor
        Vector() : m_x(constants::missing::doubleValue), m_y(constants::missing::doubleValue) {}

        /// @brief Class constructor
        ///
        /// @param[in] x The x coordinate of the vector
        /// @param[in] y The y coordinate of the vector
        Vector(const double x, const double y) : m_x(x), m_y(y) {}

        /// @brief Gets the x coordinate of the vector
        /// @returns x The x coordinate of the vector
        [[nodiscard]] double x() const
        {
            return m_x;
        }

        /// @brief Gets the x coordinate of the vector
        /// @returns x The x coordinate of the vector
        [[nodiscard]] double& x()
        {
            return m_x;
        }

        /// @brief Gets the y coordinate of the vector
        /// @returns x The y coordinate of the vector
        [[nodiscard]] double y() const
        {
            return m_y;
        }
        /// @brief Gets the y coordinate of the vector
        /// @returns x The y coordinate of the vector
        [[nodiscard]] double& y()
        {
            return m_y;
        }

        /// @brief Inplace add vector to vector.
        Vector& operator+=(const Vector& vec);

        /// @brief Inplace subtract vector from vector.
        Vector& operator-=(const Vector& vec);

        /// @brief Inplace divide vector by scalar.
        ///
        /// \note No check is made that the divisor is 0
        Vector& operator/=(const double alpha);

        /// @brief Inplace multiply vector by scalar.
        Vector& operator*=(const double alpha);

        /// @brief Normalise the vector in place.
        ///
        /// \note No check is made that the length is 0
        void normalise();

        /// @brief Compute the length of the vector
        ///
        /// @return \f$ l = \sqrt (x^2 + y^2) \f$
        double length() const;

        /// @brief Compute the length squared of the vector
        ///
        /// @return \f$ l = x^2 + y^2 \f$
        double lengthSquared() const;

    private:
        double m_x; ///< The x coordinate of the vector
        double m_y; ///< The y coordinate of the vector
    };

    /// @brief Return the normalised vector.
    ///
    /// \note No check is made that the length is 0
    Vector normalise(const Vector& vec);

    /// @brief Compute the dot product of two vectors.
    double dot(const Vector& v1, const Vector& v2);

    /// @brief Compute the angle between two vector.
    double angleBetween(const Vector& v1, const Vector& v2);

    /// @brief Add two vectors
    Vector operator+(const Vector& v1, const Vector& v2);

    /// @brief Subtract vector from another
    Vector operator-(const Vector& v1, const Vector& v2);

    /// @brief Multiply vector by a scalar
    Vector operator*(const double alpha, const Vector& vec);

    /// @brief Multiply vector by a scalar
    Vector operator*(const Vector& vec, const double alpha);

    /// @brief Divide vector by a scalar
    Vector operator/(const Vector& vec, const double alpha);

} // namespace meshkernel

inline meshkernel::Vector& meshkernel::Vector::operator+=(const Vector& vec)
{
    m_x += vec.m_x;
    m_y += vec.m_y;
    return *this;
}

inline meshkernel::Vector& meshkernel::Vector::operator-=(const Vector& vec)
{
    m_x -= vec.m_x;
    m_y -= vec.m_y;
    return *this;
}

inline meshkernel::Vector& meshkernel::Vector::operator/=(const double alpha)
{
    m_x /= alpha;
    m_y /= alpha;
    return *this;
}

inline meshkernel::Vector& meshkernel::Vector::operator*=(const double alpha)
{
    m_x *= alpha;
    m_y *= alpha;
    return *this;
}

inline void meshkernel::Vector::normalise()
{
    double lengthInv = 1.0 / length();
    m_x *= lengthInv;
    m_y *= lengthInv;
}

inline double meshkernel::Vector::length() const
{
    // TODO check implementation of hypot.
    return std::hypot(m_x, m_y);
}

inline double meshkernel::Vector::lengthSquared() const
{
    return m_x * m_x + m_y * m_y;
}

inline meshkernel::Vector meshkernel::normalise(const Vector& vec)
{
    double lengthInv = 1.0 / vec.length();
    Vector normalised(vec.x() * lengthInv, vec.y() * lengthInv);
    return normalised;
}

inline double meshkernel::dot(const Vector& v1, const Vector& v2)
{
    return v1.x() * v2.x() + v1.y() * v2.y();
}

inline double meshkernel::angleBetween(const Vector& v1, const Vector& v2)
{
    return std::atan2(v1.y() * v2.x() - v1.x() * v2.y(), v1.x() * v2.x() + v1.y() * v2.y());
}

inline meshkernel::Vector meshkernel::operator+(const Vector& v1, const Vector& v2)
{
    return Vector(v1.x() + v2.x(), v1.y() + v2.y());
}

inline meshkernel::Vector meshkernel::operator-(const Vector& v1, const Vector& v2)
{
    return Vector(v1.x() - v2.x(), v1.y() - v2.y());
}

inline meshkernel::Vector meshkernel::operator*(const double alpha, const Vector& vec)
{
    return Vector(alpha * vec.x(), alpha * vec.y());
}

inline meshkernel::Vector meshkernel::operator*(const Vector& vec, const double alpha)
{
    return alpha * vec;
}

inline meshkernel::Vector meshkernel::operator/(const Vector& vec, const double alpha)
{
    return Vector(vec.x() / alpha, vec.y() / alpha);
}
