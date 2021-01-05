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

namespace meshkernelapi
{
    /// @brief A struct used to describe the values of an unstructured mesh in a C-compatible manner
    ///
    /// \see MeshGeometryDimensions
    struct MeshGeometry
    {
        /// @brief The nodes composing each mesh 2d edge
        int* edge_nodes = nullptr;

        /// @brief The nodes composing each mesh 2d face
        int* face_nodes = nullptr;

        /// @brief The edges composing each mesh 2d face
        int* edge_faces = nullptr;

        /// @brief The left and right face indices that an edge separates
        int* face_edges = nullptr;

        /// @brief The x-coordinates of the mesh nodes
        double* nodex = nullptr;

        /// @brief The y-coordinates of the mesh nodes
        double* nodey = nullptr;

        /// @brief The z-coordinates of the mesh nodes
        double* nodez = nullptr;

        /// @brief The x-coordinates of the mesh edges middle points
        double* edgex = nullptr;

        /// @brief The y-coordinates of the mesh edges middle points
        double* edgey = nullptr;

        /// @brief  The z-coordinates of the mesh edges middle points
        double* edgez = nullptr;

        /// @brief The x-coordinates of the mesh faces mass centers
        double* facex = nullptr;

        /// @brief The y-coordinates of the mesh faces mass centers
        double* facey = nullptr;

        /// @brief The x-coordinates of the mesh faces mass centers
        double* facez = nullptr;

        /// @brief The number of mesh nodes
        int numnode;

        /// @brief The number of edges
        int numedge;

        /// @brief The number of faces
        int numface;

        /// @brief The maximum amount of nodes a face can have
        int maxnumfacenodes;
    };

} // namespace meshkernelapi
