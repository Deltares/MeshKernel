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

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Parameters.hpp>

#include <memory>

namespace meshkernel
{
    class Polygons;
    class CurvilinearGrid;

    /// @brief A class implementing the generation of a uniform curvilinear grid
    class CurvilinearGridCreateUniform
    {
    public:
        /// @brief Class constructor
        ///
        /// @param[in] MakeGridParameters The structure containing the make grid parameters
        /// @param[in] projection The projection to use
        CurvilinearGridCreateUniform(const MakeGridParameters& MakeMeshParameters, Projection projection);

        /// @brief Compute an uniform curvilinear grid using the make mesh parameters
        CurvilinearGrid Compute() const;

        /// @brief Compute an uniform curvilinear grid in one polygon
        /// @param[in] polygons The input polygons
        /// @param[in] polygonIndex The polygon index
        CurvilinearGrid Compute(std::shared_ptr<Polygons> polygons, size_t polygonIndex) const;

    private:
        static std::vector<std::vector<Point>> computeSpherical(double OriginXCoordinate,
                                                                double OriginYCoordinate,
                                                                size_t numM,
                                                                size_t numN,
                                                                double XGridBlockSize,
                                                                double YGridBlockSize,
                                                                double angle);

        static std::vector<std::vector<Point>> computeCartesian(double OriginXCoordinate,
                                                                double OriginYCoordinate,
                                                                size_t numM,
                                                                size_t numN,
                                                                double XGridBlockSize,
                                                                double YGridBlockSize,
                                                                double angle);

        Projection m_projection; ///< The projection to use

        // regular grid
        size_t m_numM;
        size_t m_numN;
        double m_angle;
        double m_XGridBlockSize;
        double m_YGridBlockSize;
        double m_OriginXCoordinate;
        double m_OriginYCoordinate;
        double m_block_size_x;
        double m_block_size_y;
    };
} // namespace meshkernel
