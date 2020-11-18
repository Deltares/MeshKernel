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
    struct MeshGeometry
    {
        int* edge_nodes = nullptr;
        int* face_nodes = nullptr;
        int* edge_faces = nullptr;
        int* face_edges = nullptr;
        int* face_links = nullptr;

        double* nnodex = nullptr;
        double* nnodey = nullptr;
        int* nedge_nodes = nullptr;
        double* nbranchlengths = nullptr;
        int* nbranchgeometrynodes = nullptr;
        double* ngeopointx = nullptr;
        double* ngeopointy = nullptr;
        int* nbranchorder = nullptr;
        int* branchidx = nullptr;
        double* branchoffsets = nullptr;

        double* nodex = nullptr;
        double* nodey = nullptr;
        double* nodez = nullptr;
        double* edgex = nullptr;
        double* edgey = nullptr;
        double* edgez = nullptr;
        double* facex = nullptr;
        double* facey = nullptr;
        double* facez = nullptr;

        double* layer_zs = nullptr;
        double* interface_zs = nullptr;
        int startIndex = 0;
    };

} // namespace meshkernelapi
