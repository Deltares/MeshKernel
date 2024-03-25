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
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>
#include <MeshKernel/UndoActions/UndoActionStack.hpp>

namespace meshkernelapi
{
    /// @brief The class holding the state of the C API library
    struct MeshKernelState
    {
        /// @brief Default constructor
        MeshKernelState() = default;

        /// @brief Constructor initializing mesh and contacts classes
        /// @param[in] projection The projection to use
        MeshKernelState(meshkernel::Projection projection) : m_projection(projection)
        {
            m_mesh1d = std::make_unique<meshkernel::Mesh1D>(projection);
            m_mesh2d = std::make_unique<meshkernel::Mesh2D>(projection);
            m_network1d = std::make_unique<meshkernel::Network1D>(projection);
            m_contacts = std::make_unique<meshkernel::Contacts>(*m_mesh1d, *m_mesh2d);
            m_curvilinearGrid = std::make_unique<meshkernel::CurvilinearGrid>(projection);
        }

        // Geometrical entities instances
        std::unique_ptr<meshkernel::Mesh1D> m_mesh1d;                   ///< Unique pointer to meshkernel::Mesh1D instance
        std::unique_ptr<meshkernel::Network1D> m_network1d;             ///< Unique pointer to meshkernel::Network1D instance
        std::unique_ptr<meshkernel::Mesh2D> m_mesh2d;                   ///< Unique pointer to meshkernel::Mesh2D instance
        std::unique_ptr<meshkernel::Contacts> m_contacts;               ///< Unique pointer to meshkernel::Contacts instance
        std::unique_ptr<meshkernel::CurvilinearGrid> m_curvilinearGrid; ///< Unique pointer to meshkernel::CurvilinearGrid instance

        // Algorithms instances (interactivity)
        std::unique_ptr<meshkernel::OrthogonalizationAndSmoothing> m_meshOrthogonalization;               ///< Unique pointer to meshkernel::OrthogonalizationAndSmoothing instance
        std::unique_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines;             ///< Unique pointer to meshkernel::CurvilinearGridFromSplines instance
        std::unique_ptr<meshkernel::CurvilinearGridOrthogonalization> m_curvilinearGridOrthogonalization; ///< Unique pointer to meshkernel::CurvilinearGridOrthogonalization instance
        std::unique_ptr<meshkernel::CurvilinearGridLineShift> m_curvilinearGridLineShift;                 ///< Unique pointer to meshkernel::CurvilinearGridLineShift instance

        // Exclusively owned state
        meshkernel::Projection m_projection{meshkernel::Projection::cartesian}; ///< Projection used by the meshes

        meshkernel::UndoActionStack m_undoStack; ///< Stack of undo actions
    };

} // namespace meshkernelapi
