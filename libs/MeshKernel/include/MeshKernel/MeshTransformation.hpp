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
#include <concepts>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Exceptions.hpp"
#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{

    /// @brief Ensure any instantiation of the MeshTransformation Compute function is with the correct operation
    template <typename Function>
    concept HasTransformationFunction = requires(Function op, Point p) {{ op(p)} -> std::same_as<Point>; };

    /// @brief Ensure any instantiation of the MeshTransformation Compute function is able to determine the projection.
    template <typename Function>
    concept HasTransformationProjection = requires(Function op) {{ op.TransformationProjection()} -> std::same_as<Projection>; };

    /// @brief To ensure the MeshTransformation Compute template parameter has a valid interface.
    template <typename Function>
    concept TransformationFunction = HasTransformationFunction<Function> && HasTransformationProjection<Function>;

    /// @brief Apply a translation transformation to a point or a vector.
    class Translation
    {
    public:
        /// @brief Default constructor, default is no translation
        Translation() = default;

        /// @brief Construct with user defined translation
        explicit Translation(const Vector& trans) : m_translation(trans) {}

        /// @brief Get the projection required for the translation
        [[nodiscard]] Projection TransformationProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief Reset translation to identity translation (i.e. no translation)
        void identity()
        {
            m_translation = {0.0, 0.0};
        }

        /// @brief Reset the translation to a new translation quantity
        void reset(const Vector& trans)
        {
            m_translation = trans;
        }

        /// @brief Get the current defined translation vector.
        const Vector& vector() const
        {
            return m_translation;
        }

        /// @brief Compose two translation objects.
        Translation compose(const Translation& trans) const
        {
            return Translation(m_translation + trans.m_translation);
        }

        /// @brief Apply the translation to a point in Cartesian coordinate system
        Point operator()(const Point& pnt) const
        {
            return pnt + m_translation;
        }

        /// @brief Apply the translation to a vector in Cartesian coordinate system
        Vector operator()(const Vector& vec) const
        {
            return vec + m_translation;
        }

    private:
        /// @brief The translation values
        Vector m_translation{0.0, 0.0};
    };

    /// @brief Apply a rotation transformation to a point or a vector.
    class Rotation
    {
    public:
        /// @brief Default constructor, default is theta = 0.
        Rotation() = default;

        /// @brief Construct with user defined rotation angle, in radians.
        explicit Rotation(const double angle) : m_theta(angle), m_cosTheta(std::cos(angle)), m_sinTheta(std::sin(angle)) {}

        /// @brief Get the projection required for the rotation
        [[nodiscard]] Projection TransformationProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief Reset rotation to identity translation (i.e. no rotation, theta = 0)
        void identity()
        {
            m_theta = 0.0;
            m_cosTheta = 1.0;
            m_sinTheta = 0.0;
        }

        /// @brief Reset the rotation to a new rotation angle
        void reset(const double angle)
        {
            m_theta = angle;
            m_cosTheta = std::cos(m_theta);
            m_sinTheta = std::sin(m_theta);
        }

        /// @brief Get the current defined rotation angle
        double angle() const
        {
            return m_theta;
        }

        /// @brief Compose two rotation objects.
        Rotation compose(const Rotation& rot) const
        {
            return Rotation(m_theta + rot.m_theta);
        }

        /// @brief Apply the rotation to a point in Cartesian coordinate system
        Point operator()(const Point& pnt) const
        {
            Point result({m_cosTheta * pnt.x - m_sinTheta * pnt.y,
                          m_sinTheta * pnt.x + m_cosTheta * pnt.y});
            return result;
        }

        /// @brief Apply the rotation to a vector in Cartesian coordinate system
        Vector operator()(const Vector& vec) const
        {
            Vector result({m_cosTheta * vec.x() - m_sinTheta * vec.y(),
                           m_sinTheta * vec.x() + m_cosTheta * vec.y()});
            return result;
        }

    private:
        /// @brief The rotation angle, theta
        double m_theta = 0.0;

        /// @brief cosine of the rotation angle
        double m_cosTheta = 1.0;

        /// @brief sine of the rotation angle
        double m_sinTheta = 0.0;
    };

    /// @brief A composition of translation and rotation transformations
    class RigidBodyTransformation
    {
    public:
        /// @brief Default constructor, default is no transformation
        RigidBodyTransformation() = default;

        /// @brief Get the projection required for the transformation
        [[nodiscard]] Projection TransformationProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief Reset transformation to identity transformation (i.e. no transformation)
        void identity()
        {
            m_rotation.identity();
            m_translation.identity();
        }

        /// @brief Compose rotation and transformation object.
        ///
        /// Will be applied: rot (transformation).
        void compose(const Rotation& rot)
        {
            m_rotation = rot.compose(m_rotation);
            m_translation.reset(rot(m_translation.vector()));
        }

        /// @brief Compose translation and transformation object.
        ///
        /// Will be applied: translation (transformation).
        void compose(const Translation& trans)
        {
            m_translation = trans.compose(m_translation);
        }

        /// @brief Get the current rotation.
        const Rotation& rotation() const
        {
            return m_rotation;
        }

        /// @brief Get the current translation.
        const Translation& translation() const
        {
            return m_translation;
        }

        /// @brief Apply the transformation to a point in Cartesian coordinate system
        Point operator()(const Point& pnt) const
        {
            Point result = m_rotation(pnt);
            result = m_translation(result);
            return result;
        }

        /// @brief Apply the transformation to a vector in Cartesian coordinate system
        Vector operator()(const Vector& vec) const
        {
            Vector result = m_rotation(vec);
            result = m_translation(result);
            return result;
        }

    private:
        /// @brief The rotation part of the transformation
        Rotation m_rotation;

        /// @brief The translation part of the transformation
        Translation m_translation;
    };

    /// @brief Apply a transformation to a mesh
    class MeshTransformation
    {
    public:
        /// @brief Apply a transformation to a mesh with a Cartesian projection
        template <TransformationFunction Transformation>
        static void Compute(Mesh& mesh, Transformation transformation)
        {
            if (mesh.m_projection != transformation.TransformationProjection())
            {
                throw MeshKernelError("Incorrect mesh coordinate system, expecting '{}', found '{}'",
                                      ToString(transformation.TransformationProjection()), ToString(mesh.m_projection));
            }

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(mesh.GetNumNodes()); ++i)
            {
                if (mesh.m_nodes[i].IsValid())
                {
                    mesh.m_nodes[i] = transformation(mesh.m_nodes[i]);
                }
            }
        }
    };

} // namespace meshkernel
