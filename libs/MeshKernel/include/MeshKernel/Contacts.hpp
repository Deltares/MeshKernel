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

/// \namespace meshkernel
/// @brief Contains the logic of the C++ static library
namespace meshkernel
{
    /// @brief A class describing an 1d-2d contacts.
    ///
    /// The responsibility of the Contacts class is connecting a 1d mesh to a 2d mesh.
    /// The class has a reference to the Mesh1D and the Mesh2D instances that will be connected.
    /// A contact is defined by the indices of the connected 1d node and 2d face.
    class Contacts
    {
    public:
        /// @brief Default constructor
        Contacts() = default;

        /// @brief Constructor taking the 1d and 2d meshes to connect
        /// @param[in] mesh1d       The mesh1d to connect
        /// @param[in] mesh2d       The mesh2d to connect
        Contacts(std::shared_ptr<Mesh1D> mesh1d,
                 std::shared_ptr<Mesh2D> mesh2d);

        /// @brief Computes 1d-2d contacts,
        /// where every single 1d node is connected to one 2d face circumcenter (ggeo_make1D2Dinternalnetlinks_dll).
        ///
        /// Each non-boundary 1d node is connected to single 2d face.
        /// The figure below shows two 2d meshes, a 1d mesh between them, and the 1d-2d contacts (in red).
        /// For the 1d nodes not overlapping a 2d mesh,
        /// a ray starting from the current node n is computed (dashed blue ray).
        /// This ray is normal to the segment connecting the previous (n-1) and next one 1d node (n+1),
        /// the connecting segment is shown with a green dashed line.
        /// The ray is extended for 5 times the length of the connecting segment.
        /// The current 1d node is connected to the first boundary 2d face crossing the ray,
        /// first in the left direction and then in the right direction.
        /// By doing so a 1d mesh can be connected on the left and right sides of a mesh 2d boundary,
        /// for example when the 1d part represents a river and the 2d part the river banks.
        /// The 1d nodes overlapping the 2d mesh are directly connected to the face including them.
        /// \image html ComputeSingleContacts.jpg  "1d mesh connecting to 2d mesh using the ComputeSingleContacts algorithm. Contacts are shown in red."
        ///
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)
        /// @param[in] polygons     The polygons selecting the area where the 1d-2d contacts will be generated.
        /// @param[in] projectionFactor     The projection factor used for generating contacts when 1d nodes are not inside mesh2d
        void ComputeSingleContacts(const std::vector<bool>& oneDNodeMask, const Polygons& polygons, double projectionFactor);

        /// @brief Computes 1d-2d contacts,
        /// where a single 1d node is connected to multiple 2d face circumcenters (ggeo_make1D2Dembeddedlinks_dll).
        ///
        /// Each internal 1d node is connected to multiple 2d faces.
        /// This type of contacts should be used when the lengths of the 1d mesh edges are considerably larger
        /// than the 2d mesh edges and generating a single contact for each 1d node is not sufficient.
        /// In this algorithm, only the internal 1d nodes are connected.
        /// The following figure shows a 1d mesh overlapping a 2d mesh.
        /// For the node n, the closest 2d faces within a search radius are found
        /// and it is determined if those faces cross are crossed by the current 1d edge starting at node n and ending at and n+1.
        /// If the answer is positive, a contact is generated between the face
        /// and the closest 1d node composing the current 1d edge (i.e. n or n+1).
        /// The procedure is repeated for each 1d node.
        /// \image html ComputeMultipleContacts.jpg  "1d mesh connecting to 2d mesh using the ComputeMultipleContacts algorithm. Contacts are shown in red."
        ///
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)
        void ComputeMultipleContacts(const std::vector<bool>& oneDNodeMask);

        /// @brief Computes 1d-2d contacts,
        /// where a 2d face per polygon is connected to the closest 1d node (ggeo_make1D2Droofgutterpipes_dll).
        ///
        /// The algorithms works as follows:
        /// - Find the 2d face within each polygon closest to a 1d node.
        /// - Per polygon create one contact from the 2d circumcenters to the 1d node.
        /// \image html ComputeContactsWithPolygons.svg  "1d mesh connecting to 2d mesh using the ComputeContactsWithPolygons algorithm. Contacts are shown in red. Polygons in green."
        ///
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)
        /// @param[in] polygons     The polygons to connect (Polygons class can have multiple polygons)
        void ComputeContactsWithPolygons(const std::vector<bool>& oneDNodeMask, const Polygons& polygons);

        /// @brief Computes 1d-2d contacts,
        /// where 1d nodes are connected to the 2d faces mass centers containing the input point (ggeo_make1D2Dstreetinletpipes_dll).
        ///
        /// With this algorithm, each 2d face containing a point is connected to the 1d node closest to point itself.
        /// The search of the 2d faces and the closest 1d nodes uses RTrees.
        /// \image html ComputeContactsWithPoints.jpg  "2d faces containing the input points connecting to the 1d mesh. Contacts are shown in red and the input points in blue."
        ///
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)
        /// @param[in] points       The points selecting the faces to connect
        void ComputeContactsWithPoints(const std::vector<bool>& oneDNodeMask, const std::vector<Point>& points);

        /// @brief Computes 1d-2d contacts,
        /// where 1d nodes are connected to the closest 2d faces at the boundary (ggeo_make1D2DRiverLinks_dll).
        ///
        /// The algorithms works as follows:
        /// - For each 1d node, find the closest 2d boundary faces within the search radius.
        /// - If a boundary face can be connected to multiple oned nodes, choose the closest one.
        /// - Generate the 1d-2d contacts.Index m_numM = 0;                                     ///< Number of columns in the curvilinear grid
        /// \image html ComputeBoundaryContacts.jpg  "1d mesh connecting to 2d mesh using the ComputeBoundaryContacts algorithm. Contacts are shown in red.
        /// The mesh 2d boundary faces are connected to the closest 1d nodes."
        ///
        /// @param[in] oneDNodeMask The mask to apply to 1d nodes (true = connect node, false = do not generate contacts)
        /// @param[in] polygons     The polygons selecting the area where the 1d-2d contacts will be generated.
        /// @param[in] searchRadius The radius used for searching neighboring faces, if equal to constants::missing::doubleValue, the search radius will be calculated internally.
        void ComputeBoundaryContacts(const std::vector<bool>& oneDNodeMask,
                                     const Polygons& polygons,
                                     double searchRadius);

        /// @brief Gets the 1d mesh indices
        /// @return Vector of 1d mesh indices
        std::vector<UInt> const& Mesh1dIndices() const { return m_mesh1dIndices; }

        /// @brief Gets the 2d mesh indices
        /// @return Vector of 2d mesh indices
        std::vector<UInt> const& Mesh2dIndices() const { return m_mesh2dIndices; }

        /// @brief Sets the 1d and 2d mesh indices
        /// @param[in] mesh1dIndices The 1d mesh indices
        /// @param[in] mesh2dIndices The 2d mesh indices
        void SetIndices(const std::vector<meshkernel::UInt>& mesh1dIndices,
                        const std::vector<meshkernel::UInt>& mesh2dIndices)
        {
            m_mesh1dIndices = mesh1dIndices;
            m_mesh2dIndices = mesh2dIndices;
            m_areComputed = true;
        }

        /// @brief checks whether contacts have been computed
        /// @return True if the contact is crossing an existing contact
        [[nodiscard]] bool AreComputed() const { return m_areComputed; };

    private:
        /// @brief Asserts if a contact is crossing a 1d mesh edge
        /// @param[in] node The 1d node index (start of the contact)
        /// @param[in] face The 2d face index (end of the contact)
        /// @return True if the contact is crossing a 1d mesh edge
        [[nodiscard]] bool IsContactIntersectingMesh1d(UInt node, UInt face) const;

        /// @brief Asserts if a contact is crossing an existing contact
        /// @param[in] node The 1d node index (start of the contact)
        /// @param[in] face The 2d face index (end of the contact)
        /// @return True if the contact is crossing an existing contact
        [[nodiscard]] bool IsContactIntersectingContact(UInt node, UInt face) const;

        /// @brief Connect the current 1D line segment with the faces that intersect a semiline originating from the current node and perpendicular to the current 1D edge.
        /// @param[in] node The 1d node index (start of the contact)
        /// @param[in] projectionFactor The semiline length, as a multiplier of the current ad edge length
        void Connect1dNodesWithCrossingFaces(UInt node,
                                             double projectionFactor);

        std::shared_ptr<Mesh1D> m_mesh1d;  ///< The 1-d mesh to connect
        std::shared_ptr<Mesh2D> m_mesh2d;  ///< The 2-d mesh to connect
        std::vector<UInt> m_mesh1dIndices; ///< The indices of the connected 1-d nodes
        std::vector<UInt> m_mesh2dIndices; ///< The indices of the connected 2-d faces
        bool m_areComputed = false;        ///< Indicates whether contacts have been computed
    };
} // namespace meshkernel
