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

#include "ApiCache/MeshBoundariesAsPolygonCache.hpp"

#include <memory>

#include "MeshKernel/Contacts.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridOrthogonalization.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridSmoothing.hpp"
#include "MeshKernel/Mesh1D.hpp"
#include "MeshKernel/Mesh2D.hpp"
#include "MeshKernel/OrthogonalizationAndSmoothing.hpp"

#include "MeshKernelApi/ApiCache/CurvilinearBoundariesAsPolygonCache.hpp"
#include "MeshKernelApi/ApiCache/FacePolygonPropertyCache.hpp"
#include "MeshKernelApi/ApiCache/HangingEdgeCache.hpp"
#include "MeshKernelApi/ApiCache/NodeInPolygonCache.hpp"
#include "MeshKernelApi/ApiCache/ObtuseTriangleCentreCache.hpp"
#include "MeshKernelApi/ApiCache/PolygonRefinementCache.hpp"
#include "MeshKernelApi/ApiCache/SmallFlowEdgeCentreCache.hpp"

namespace meshkernelapi
{

    /// @brief The class holding the state of the C API library
    struct MeshKernelState
    {

        /// @brief Default constructor
        MeshKernelState() = default;

        /// @brief Simple constructor
        explicit MeshKernelState(meshkernel::Projection projection) : m_projection(projection)
        {
            m_mesh1d = std::make_shared<meshkernel::Mesh1D>(projection);
            m_mesh2d = std::make_shared<meshkernel::Mesh2D>(projection);
            m_network1d = std::make_shared<meshkernel::Network1D>(projection);
            m_contacts = std::make_shared<meshkernel::Contacts>(*m_mesh1d, *m_mesh2d);
            m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(projection);
            m_frozenLinesCounter = 0;
        }

        // Geometrical entities instances
        std::shared_ptr<meshkernel::Mesh1D> m_mesh1d;                   ///< Shared pointer to meshkernel::Mesh1D instance
        std::shared_ptr<meshkernel::Network1D> m_network1d;             ///< Shared pointer to meshkernel::Network1D instance
        std::shared_ptr<meshkernel::Mesh2D> m_mesh2d;                   ///< Shared pointer to meshkernel::Mesh2D instance
        std::shared_ptr<meshkernel::Contacts> m_contacts;               ///< Shared pointer to meshkernel::Contacts instance
        std::shared_ptr<meshkernel::CurvilinearGrid> m_curvilinearGrid; ///< Shared pointer to meshkernel::CurvilinearGrid instance

        // Algorithms instances (interactivity)
        std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing> m_meshOrthogonalization;                  ///< Shared pointer to meshkernel::OrthogonalizationAndSmoothing instance
        std::shared_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines;                ///< Shared pointer to meshkernel::CurvilinearGridFromSplines instance
        std::shared_ptr<meshkernel::CurvilinearGridLineShift> m_curvilinearGridLineShift;                    ///< Shared pointer to meshkernel::CurvilinearGridLineShift instance
        std::unordered_map<meshkernel::UInt, std::pair<meshkernel::Point, meshkernel::Point>> m_frozenLines; ///< Map for string the frozen lines
        meshkernel::UInt m_frozenLinesCounter = 0;                                                           ///< An increasing counter for returning the id of frozen lines to the client

        // Exclusively owned state
        meshkernel::Projection m_projection{meshkernel::Projection::cartesian}; ///< Projection used by the meshes

        // Cached values, used when dimensions are computed first, followed by values being retrieved in a separate call
        std::shared_ptr<FacePolygonPropertyCache> m_facePropertyCache;                   ///< face property cache
        std::shared_ptr<CurvilinearBoundariesAsPolygonCache> m_boundariesAsPolygonCache; ///< boundaries as polygon cache
        std::shared_ptr<MeshBoundariesAsPolygonCache> m_meshBoundariesAsPolygonCache;    ///< boundaries as polygon cache
        std::shared_ptr<PolygonRefinementCache> m_polygonRefinementCache;                ///< polygon refinement cache
        std::shared_ptr<NodeInPolygonCache> m_nodeInPolygonCache;                        ///< node in polygon cache
        std::shared_ptr<SmallFlowEdgeCentreCache> m_smallFlowEdgeCentreCache;            ///< small flow edge centres cache
        std::shared_ptr<HangingEdgeCache> m_hangingEdgeCache;                            ///< hanging edge id cache
        std::shared_ptr<ObtuseTriangleCentreCache> m_obtuseTriangleCentreCache;          ///< centre of obtuse triangles cache
    };

} // namespace meshkernelapi
