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

namespace meshkernelapi
{
    /// @brief A struct used to describe the geometry dimensions of an unstructured mesh in a C-compatible manner
    ///
    /// \see MeshGeometry
    struct MeshGeometryDimensions
    {
        /// @brief name
        char name[255];

        /// @brief The mesh dimension (e.g. 1D or 2D)
        int dim;

        /// @brief The number of mesh nodes
        int numnode;

        /// @brief The number of edges
        int numedge;

        /// @brief The number of faces
        int numface;

        /// @brief The maximum amount of nodes a face can have
        int maxnumfacenodes;

        /// @brief The index of the layer the current mesh represents
        int numlayer;

        /// @brief The layer type the current mesh represents
        int layertype;

        /// @brief The number of network nodes
        int nnodes;

        /// @brief The number of network branches
        int nbranches;

        /// @brief The number of geometry nodes
        int ngeometry;

        /// @brief The epsg code of the mesh projection
        int epsg;
    };

} // namespace meshkernelapi
