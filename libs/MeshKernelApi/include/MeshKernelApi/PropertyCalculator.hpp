//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include "MeshKernel/SampleInterpolator.hpp"

#include "MeshKernelApi/GeometryList.hpp"
#include "MeshKernelApi/State.hpp"

namespace meshkernelapi
{

    /// @brief Base class for calculating properties for a mesh
    class PropertyCalculator
    {
    public:
        /// @brief Destructor
        virtual ~PropertyCalculator() = default;

        /// @brief Determine is the calculator can compute the desired results correctly.
        ///
        /// This has a default of checking that the mesh2d is not null and the number of nodes is greater than zero.
        virtual bool IsValid(const MeshKernelState& state) const;

        /// @brief Calculate the property
        virtual void Calculate(const MeshKernelState& state, const GeometryList& geometryList) const = 0;

        /// @brief Determine the size of the vector required to store the calculated properties
        virtual int Size(const MeshKernelState& state) const = 0;
    };

    /// @brief Calculator for orthogonality of a mesh.
    class OrthogonalityPropertyCalculator : public PropertyCalculator
    {
    public:
        /// @brief Calculate the orthogonality for a mesh
        void Calculate(const MeshKernelState& state, const GeometryList& geometryList) const override;

        /// @brief Determine the size of the orthogonality vector required
        int Size(const MeshKernelState& state) const override;
    };

    /// @brief Calculator for the edge lengths for a mesh
    class EdgeLengthPropertyCalculator : public PropertyCalculator
    {
    public:
        /// @brief Calculate the edge-length for a mesh
        void Calculate(const MeshKernelState& state, const GeometryList& geometryList) const override;

        /// @brief Determine the size of the edge-length vector required
        int Size(const MeshKernelState& state) const override;
    };

    /// @brief Interpolate the depths at the mesh node points.
    class InterpolatedSamplePropertyCalculator : public PropertyCalculator
    {
    public:
        /// @brief Constructor
        InterpolatedSamplePropertyCalculator(const GeometryList& sampleData,
                                             const meshkernel::Projection projection,
                                             const int propertyId);

        /// @brief Determine is the calculator can interpolate depth values correctly
        bool IsValid(const MeshKernelState& state) const override;

        /// @brief Calculate the edge-length for a mesh
        void Calculate(const MeshKernelState& state, const GeometryList& geometryList) const override;

        /// @brief Determine the size of the edge-length vector required
        int Size(const MeshKernelState& state) const override;

    private:
        /// @brief Interpolator for the samples
        std::unique_ptr<meshkernel::SampleInterpolator> m_sampleInterpolator;

        /// @brief Projection sued for sample data.
        meshkernel::Projection m_projection;

        /// @brief Property id.
        int m_propertyId = -1;
    };

} // namespace meshkernelapi
