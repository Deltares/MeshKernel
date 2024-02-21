//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <memory>

#include "MeshKernel/BaseMeshUndoAction.hpp"
#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Point.hpp"

namespace meshkernel
{
    class Mesh;

    class AddEdgeAction : public BaseMeshUndoAction<AddEdgeAction, Mesh>
    {
    public:
        static std::unique_ptr<AddEdgeAction> Create(Mesh& mesh, const UInt id, const UInt start, const UInt end);

        AddEdgeAction(Mesh& mesh, const UInt id, const UInt start, const UInt end);

        UInt EdgeId() const;

        const Edge& GetEdge() const;

    private:
        UInt m_edgeId;
        Edge m_edge;
    };

} // namespace meshkernel
