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

#include <memory>

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>

namespace meshkernelapi
{
    /// @brief The class holding the state of the C API library
    class MeshKernelState
    {
    public:
        /// @brief Default constructor
        /// @returns
        MeshKernelState() = default;

        /// @brief Constructor initializing mesh and contacts classes
        /// @param[in] projection The projection to use
        MeshKernelState(meshkernel::Projection projection) : m_projection(projection)
        {
            m_mesh1d = std::make_shared<meshkernel::Mesh1D>();
            m_mesh2d = std::make_shared<meshkernel::Mesh2D>();
            m_contacts = std::make_shared<meshkernel::Contacts>(m_mesh1d, m_mesh2d);
            m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>();
        }

        // Geometrical entities instances
        std::shared_ptr<meshkernel::Mesh1D> m_mesh1d;                   ///< Shared pointer to meshkernel::Mesh1D instance
        std::shared_ptr<meshkernel::Mesh2D> m_mesh2d;                   ///< Shared pointer to meshkernel::Mesh2D instance
        std::shared_ptr<meshkernel::Contacts> m_contacts;               ///< Shared pointer to meshkernel::Contacts instance
        std::shared_ptr<meshkernel::CurvilinearGrid> m_curvilinearGrid; ///< Shared pointer to meshkernel::CurvilinearGrid instance

        // Algorithms instances (interactivity)
        std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing> m_meshOrthogonalization;               ///< Shared pointer to meshkernel::OrthogonalizationAndSmoothing instance
        std::shared_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines;             ///< Shared pointer to meshkernel::CurvilinearGridFromSplines instance
        std::shared_ptr<meshkernel::CurvilinearGridOrthogonalization> m_curvilinearGridOrthogonalization; ///< Shared pointer to meshkernel::CurvilinearGridOrthogonalization instance
        std::shared_ptr<meshkernel::CurvilinearGridLineShift> m_curvilinearGridLineShift;                 ///< Shared pointer to meshkernel::CurvilinearGridLineShift instance

        // Exclusively owned state
        meshkernel::Projection m_projection; ///< Projection used by the meshes
    };

} // namespace meshkernelapi
