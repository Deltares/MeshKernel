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
        /// @brief Default constructor
        Contacts() = default;

        /// @brief Constructor taking the 1d and 2d meshes to connect
        /// @param[in] mesh1d The mesh1d to connect
        /// @param[in] mesh2d The mesh2d to connect
        /// @param[in] oneDNodeMask The mask used for masking 1d nodes
        Contacts(std::shared_ptr<Mesh1D> mesh1d, std::shared_ptr<Mesh2D> mesh2d, const std::vector<bool>& oneDNodeMask);

        /// @brief Computes 1d-2d connections, where every single 1d node is connected to one 2d face circumcenter (ggeo_make1D2Dinternalnetlinks_dll)
        /// @param[in] polygons The polygons where the 1d-2d connections are generated
        void ComputeSingleConnections(const Polygons& polygons);

        /// @brief Computes 1d-2d connections, where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll)
        void ComputeMultipleConnections();

        /// @brief Computes 1d-2d connections, where a 1d node is connected to the closest polygon (ggeo_make1D2Droofgutterpipes_dll)
        /// @param[in] polygons The polygons to connect (Polygons class can have multiple polygons)
        void ComputeConnectionsWithPolygons(const Polygons& polygons);

        /// @brief Computes 1d-2d connections, where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeo_make1D2Dstreetinletpipes_dll)
        /// @param[in] points the points to connect
        void ComputeConnectionsWithPoints(const std::vector<Point>& points);

        /// @brief Computes 1d-2d connections, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        void ComputeBoundaryConnections();

    private:
        /// @brief Asserts if a connection is crossing a 1d mesh edge
        /// @param node[in] The 1d node index (start of the connection)
        /// @param face[in] The 2d face index (end of the connection)
        /// @return True if the connection is crossing a 1d mesh edge
        bool IsConnectionIntersectingMesh1d(size_t node, size_t face) const;

        /// @brief Asserts if a connection is crossing an existing connection
        /// @param node[in] The 1d node index (start of the connection)
        /// @param face[in] The 2d face index (end of the connection)
        /// @return True if the connection is crossing an existing connection
        bool IsContactIntersectingContact(size_t node, size_t face) const;

        std::shared_ptr<Mesh2D> m_mesh2d;
        std::shared_ptr<Mesh1D> m_mesh1d;
        // nodes
        std::vector<size_t> m_mesh2dIndices;
        std::vector<size_t> m_mesh1dIndices;
        std::vector<bool> m_oneDNodeMask;
    };
} // namespace meshkernel