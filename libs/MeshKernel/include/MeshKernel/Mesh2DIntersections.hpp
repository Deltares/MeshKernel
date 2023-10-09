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
#include <vector>

#include "MeshKernel/Definitions.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Polygons.hpp"

namespace meshkernel
{
    /// @brief Compute the intersections of polygon inner and outer perimeters
    ///
    /// @note Uses a breadth first algorithm to reduce runtime complexity
    class Mesh2DIntersections final
    {
    public:
        /// @brief compute the edges and faces intersected by a polygon, with additional information on the intersections
        void Compute(Mesh2D& mesh, const Polygons& polygon);

        const auto& EdgeIntersections() const { return m_edgesIntersections; }
        const auto& FaceIntersections() const { return m_faceIntersections; }

    private:
        /// @brief Gets the intersection from a single polyline
        /// @param[in] polyLine An input polyline
        void GetPolylineIntersection(const Mesh2D& mesh, const std::vector<Point>& polyLine);

        std::tuple<UInt, UInt> GetIntersectionSeed(const Mesh2D& mesh, const std::vector<Point>& polyLine, const std::vector<bool>& vistedEdges) const;

        std::vector<EdgeMeshPolylineIntersection> m_edgesIntersectionsCache; ///< A cache for saving the edge intersections of one inner or outer 
        std::vector<FaceMeshPolylineIntersection> m_facesIntersectionsCache; ///< A cache for saving the local face intersections of one inner or outer 
        std::vector<EdgeMeshPolylineIntersection> m_edgesIntersections;      ///< A vector collecting all edge intersection results
        std::vector<FaceMeshPolylineIntersection> m_faceIntersections;       ///< A vector collecting all face intersection results
    };

} // namespace meshkernel
