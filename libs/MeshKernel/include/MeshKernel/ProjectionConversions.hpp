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

#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Point.hpp"
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{

    /// @brief Ensure any instantiation of the MeshTransformation Compute function is with the correct operation
    template <typename Operation>
    concept TransformationOperation = requires(Operation op, Point p)
    {
        {op(p)} -> std::same_as<Point>;
    };

    /// @brief Ensure any instantiation of the MeshConversion Compute function is able to determine source and target projection types.
    template <typename Operation>
    concept TransformationProjection = requires(Operation op)
    {
        {op.SourceProjection ()} -> std::same_as<Projection>;
        {op.TargetProjection ()} -> std::same_as<Projection>;
    };

    /// @brief Ensure the MeshConversion template parameter has a valid interface.
    template <typename Operation>
    concept ConversionFunctor = TransformationOperation<Operation> && TransformationProjection<Operation>;

    /// @brief Converts points from spherical to Cartesian coordinate system.
    class ConvertSphericalToCartesian
    {
    public :
        /// @brief Construct ConvertSphericalToCartesian with the origin of the target mesh
        ConvertSphericalToCartesian(const Point& o) : origin (o) {}

        /// @brief The coordinate system of the point parameter to the conversion operation.
        Projection SourceProjection() const
        {
            return Projection::spherical;
        }

        /// @brief The coordinate system of the point result of the conversion operation.
        Projection TargetProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief Apply the conversion of a point in Spherical coordinate system to Cartesian
        Point operator()(const Point& pnt) const
        {
            Point result = origin + Vector(GetDx(origin, pnt, Projection::spherical), GetDy(origin, pnt, Projection::spherical));
            return result;
        }

    private :
        /// @brief The origin of the target mesh.
        Point origin;
    };

    /// @brief Apply a conversion to nodes of a mesh
    class MeshConversion
    {
    public:
        /// @brief Apply a conversion to nodes of a mesh.
        template <ConversionFunctor Conversion>
        static void Compute(const Mesh& sourceMesh, Mesh& targetMesh, const Conversion& conversion)
        {
            if (sourceMesh.m_projection != conversion.SourceProjection ())
            {
                throw MeshKernelError("Incorrect source mesh coordinate system, expecting '{}', found '{}'",
                                      ToString(conversion.SourceProjection ()), ToString(sourceMesh.m_projection));
            }

            if (targetMesh.m_projection != conversion.TargetProjection ())
            {
                throw MeshKernelError("Incorrect target mesh coordinate system, expecting '{}', found '{}'",
                                      ToString(conversion.TargetProjection ()), ToString(targetMesh.m_projection));
            }

            if (sourceMesh.GetNumNodes () != targetMesh.GetNumNodes ())
            {
                throw MeshKernelError("Source and target meshes have different numbers of nodes, source = '{}', target = '{}'",
                                      sourceMesh.GetNumNodes (), targetMesh.GetNumNodes ());
            }

#pragma omp parallel for
            for (int i = 0; i < static_cast<int>(sourceMesh.GetNumNodes()); ++i)
            {
                // TODO SHould the conversion function take care of invalid nodes?
                if (sourceMesh.m_nodes[i].IsValid())
                {
                    targetMesh.m_nodes[i] = conversion(sourceMesh.m_nodes[i]);
                }
                else
                {
                    targetMesh.m_nodes[i] = sourceMesh.m_nodes[i];
                }
            }

            // TODO Does this need a parameter indicating that admininstrate needs to be called
            targetMesh.Administrate();
        }
    };

}
