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

#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/Polygons.hpp"
#include "MeshKernel/Utilities/RTreeBase.hpp"

namespace meshkernel
{
    class CasulliRefinement
    {
    public :
        static void Compute (Mesh2D& mesh, const MeshRefinementParameters& meshRefinementParameters);

    private :

        using LinkNodes = std::array<UInt, 4>;

        static void ComputeNewNodes (Mesh2D& mesh, std::vector<LinkNodes>& newNodes);

        static void StoreNewNode (const Mesh2D& mesh, const UInt knode, const UInt link1, const UInt link2, const UInt knew, std::vector<LinkNodes>& newNodes);

        static UInt IsStartEnd (const Mesh2D& mesh, const UInt nodeId, const UInt edgeId);

        static UInt IsLeftRight (const Mesh2D& mesh, const UInt elementId, const UInt edgeId);

        static UInt FindCommon (const Mesh2D& mesh, const UInt l1, const UInt l2);

    };

} // namespace meshkernel
