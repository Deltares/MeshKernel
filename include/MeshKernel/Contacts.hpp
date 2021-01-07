//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2020.
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

#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Polygons.hpp>

#include <memory>
#include <vector>

namespace meshkernel
{
    class Contacts
    {
    public:
        Contacts() = default;

        /// @brief
        /// @param mesh1d
        /// @param mesh
        /// @param projection
        Contacts(std::shared_ptr<Mesh1D> mesh1d, std::shared_ptr<Mesh2D> mesh);

        /// @brief Computes 1D-2D connections, where every single 1d node is connected to one 2d face mass center (ggeo_make1D2Dinternalnetlinks_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeSingleConnections(const meshkernel::Polygons& polygons);

        /// @brief  Computes 1D-2D connections, where a single 1d point is connected to multiple 2d face mass centers (ggeo_make1D2Dembeddedlinks_dll)
        /// @param mesh1d
        /// @param mesh
        void ComputeMultipleConnections();

        /// @brief Computes 1D-2D connections, where a 1d point is connected to the closest 2d face in polygons (ggeo_make1D2Droofgutterpipes_dll)
        /// @param mesh1d
        /// @param mesh
        /// @param polygons
        void ComputeConnectionsWithPolygons(const Polygons& polygons);

        /// @brief Computes 1D-2D connections, where 1d nodes are connected to the 2d faces mass centers containing the input points (ggeo_make1D2Dstreetinletpipes_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeConnectionsWithPoints(const std::vector<Point>& points);

        /// @brief Computes 1D-2D connections, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        /// @param mesh1d
        /// @param mesh
        /// @return
        void ComputeBoundaryConnections();

    private:
        bool IsConnectionIntersectingMesh1d(size_t node, size_t face) const;

        bool IsContactIntersectingContact(size_t node, size_t face) const;

        std::shared_ptr<Mesh2D> m_mesh2d;
        std::shared_ptr<Mesh1D> m_mesh1d;
        // nodes
        std::vector<size_t> m_mesh2dIndices;
        std::vector<size_t> m_mesh1dIndices;
    };
} // namespace meshkernel