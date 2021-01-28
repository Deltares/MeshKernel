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

#include <memory>
#include <vector>

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class describing an 1d-2d contacts
    ///
    /// The responsibility of the Contacts class is connecting a 1d mesh to a 2d mesh.
    /// The class has a reference to the Mesh1D and the Mesh2D instances that will be connected.
    /// A connection is defined by the indices of the connected 1d node and 2d face.
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
        ///
        /// Each non-boundary 1d node is connected to single 2d face.
        /// The figure below shows two 2d meshes, a 1d mesh between them, and the 1d-2d connections (in red).
        /// The boundary nodes of the 1d mesh (those sharing only one 1d edge) are not connected to any 2d face.
        /// For the 1d nodes not overlapping a 2d mesh, a ray starting from the current node n is computed (dashed blue ray).
        /// This ray is normal to the segment connecting the previous (n-1) and next one 1d node (n+1),
        /// the connecting segment is shown with a green dashed line.
        /// The ray is extended for 5 times the length of the connecting segment.
        /// The current 1d node is connected to the first boundary 2d face crossing the ray,
        /// first in the left direction and then in the right direction.
        /// By doing so a 1d mesh can be connected on the left and right sides of a mesh 2d boundary,
        /// for example when the 1d part represents a river and the 2d part the river banks.
        /// The 1d nodes overlapping the 2d mesh are directly connected to the face including them.
        /// \image html ComputeSingleConnections.jpg  "1d mesh connecting to 2d mesh using the ComputeSingleConnections algorithm. Connections are shown in red."
        /// @param[in] polygons The polygons where the 1d-2d connections are generated
        void ComputeSingleConnections(const Polygons& polygons);

        /// @brief Computes 1d-2d connections, where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll)
        ///
        /// Each internal 1d node is connected to multiple 2d faces.
        /// This type of connections should be used when the lengths of the 1d mesh edges are considerably larger
        /// than the 2d mesh edges and generating a single connection for each 1d node is not representative.
        /// In this algorithm, only the internal 1d nodes are connected.
        /// The following figure shows a 1d mesh overlapping a 2d mesh.
        /// For the node n, the closest 2d faces within a search radius are found
        /// and it is determined if those faces cross are crossed by the current 1d edge starting at node n and ending at and n+1.
        /// If the answer is positive, a connection is generated between the face
        /// and the closest 1d node composing the current 1d edge (i.e. n or n+1).
        /// The procedure is repeated for each 1d node.
        /// \image html ComputeMultipleConnections.jpg  "1d mesh connecting to 2d mesh using the ComputeMultipleConnections algorithm. Connections are shown in red."
        void ComputeMultipleConnections();

        /// @brief Computes 1d-2d connections, where a 1d node is connected to the closest polygon (ggeo_make1D2Droofgutterpipes_dll)
        ///
        /// \image html ComputeConnectionsWithPolygons.svg  "1d mesh connecting to 2d mesh using the ComputeConnectionsWithPolygonss algorithm. Connections are shown in red. Polygons in green."
        /// @param[in] polygons The polygons to connect (Polygons class can have multiple polygons)
        void ComputeConnectionsWithPolygons(const Polygons& polygons);

        /// @brief Computes 1d-2d connections, where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeso_make1D2Dstreetinletpipes_dll)
        /// @param[in] points The points to connect
        void ComputeConnectionsWithPoints(const std::vector<Point>& points);

        /// @brief Computes 1d-2d connections, where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll)
        void ComputeBoundaryConnections();

        std::vector<size_t> m_mesh2dIndices; ///< The indices of the connected 2-d faces
        std::vector<size_t> m_mesh1dIndices; ///< The indices of the connected 1-d nodes

    private:
        /// @brief Asserts if a connection is crossing a 1d mesh edge
        /// @param[in] node The 1d node index (start of the connection)
        /// @param[in] face The 2d face index (end of the connection)
        /// @return True if the connection is crossing a 1d mesh edge
        [[nodiscard]] bool IsConnectionIntersectingMesh1d(size_t node, size_t face) const;

        /// @brief Asserts if a connection is crossing an existing connection
        /// @param[in] node The 1d node index (start of the connection)
        /// @param[in] face The 2d face index (end of the connection)
        /// @return True if the connection is crossing an existing connection
        [[nodiscard]] bool IsContactIntersectingContact(size_t node, size_t face) const;

        std::shared_ptr<Mesh1D> m_mesh1d; ///< The 1-d mesh to connect
        std::shared_ptr<Mesh2D> m_mesh2d; ///< The 2-d mesh to connect
        std::vector<bool> m_oneDNodeMask; ///< The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)

        /// @brief Connect a 1d node with the face crossed by the projected normal originating from the node itself
        /// @param[in] node The 1d node index
        /// @param[in] distanceFactor The factor determining the length and the direction of the projected normal (positive right normal, negative left normal)
        void Connect1dNodesWithCrossingFace(size_t node, double distanceFactor);
    };
} // namespace meshkernel
