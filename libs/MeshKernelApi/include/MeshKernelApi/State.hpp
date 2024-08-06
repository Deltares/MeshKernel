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
#include <string>

#include <MeshKernel/Contacts.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridFromSplines.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridLineShift.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/OrthogonalizationAndSmoothing.hpp>

namespace meshkernelapi
{

    /// @brief The class holding the state of the C API library
    struct MeshKernelState
    {

        /// @brief Indicator of the current state of the MeshKernelState
        enum class CurrentState
        {
            Uninitialised, ///< MeshKernelState is currently not initialised.
            ValidMesh,     ///< MeshKernelState currently points to a valid mesh, this may be mesh1d, mesh2d or curvilinear grid.
            DeletedMesh    ///< MeshKernelState currently points to a deleted mesh.
        };

        /// @brief Default constructor
        MeshKernelState() = default;

        // MeshKernelState(MeshKernelState& mk) : m_mesh1d(mk.m_mesh1d),
        //                                        m_network1d(mk.m_network1d),
        //                                        m_mesh2d(mk.m_mesh2d),
        //                                        m_contacts(mk.m_contacts),
        //                                        m_curvilinearGrid(mk.m_curvilinearGrid),
        //                                        m_meshOrthogonalization(mk.m_meshOrthogonalization),
        //                                        m_curvilinearGridFromSplines(mk.m_curvilinearGridFromSplines),
        //                                        m_curvilinearGridOrthogonalization(mk.m_curvilinearGridOrthogonalization),
        //                                        m_curvilinearGridLineShift(mk.m_curvilinearGridLineShift),
        //                                        m_projection(mk.m_projection),
        //                                        m_state(mk.m_state)

        // {
        // }

        // MeshKernelState(MeshKernelState&& mk) : m_mesh1d(std::move(mk.m_mesh1d)),
        //                                         m_network1d(std::move(mk.m_network1d)),
        //                                         m_mesh2d(std::move(mk.m_mesh2d)),
        //                                         m_contacts(std::move(mk.m_contacts)),
        //                                         m_curvilinearGrid(std::move(mk.m_curvilinearGrid)),
        //                                         m_meshOrthogonalization(std::move(mk.m_meshOrthogonalization)),
        //                                         m_curvilinearGridFromSplines(std::move(mk.m_curvilinearGridFromSplines)),
        //                                         m_curvilinearGridOrthogonalization(std::move(mk.m_curvilinearGridOrthogonalization)),
        //                                         m_curvilinearGridLineShift(std::move(mk.m_curvilinearGridLineShift)),
        //                                         m_projection(mk.m_projection),
        //                                         m_state(mk.m_state)

        // {
        // }

        // /// @brief Constructor initializing mesh and contacts classes
        // /// @param[in] projection The projection to use
        // MeshKernelState(meshkernel::Projection projection) : m_projection(projection)
        // {
        //     m_mesh1d = std::make_unique<meshkernel::Mesh1D>(projection);
        //     m_mesh2d = std::make_unique<meshkernel::Mesh2D>(projection);
        //     m_network1d = std::make_unique<meshkernel::Network1D>(projection);
        //     m_contacts = std::make_unique<meshkernel::Contacts>(*m_mesh1d, *m_mesh2d);
        //     m_curvilinearGrid = std::make_unique<meshkernel::CurvilinearGrid>(projection);
        // }

        // // Geometrical entities instances
        // std::unique_ptr<meshkernel::Mesh1D> m_mesh1d;                   ///< Unique pointer to meshkernel::Mesh1D instance
        // std::unique_ptr<meshkernel::Network1D> m_network1d;             ///< Unique pointer to meshkernel::Network1D instance
        // std::unique_ptr<meshkernel::Mesh2D> m_mesh2d;                   ///< Unique pointer to meshkernel::Mesh2D instance
        // std::unique_ptr<meshkernel::Contacts> m_contacts;               ///< Unique pointer to meshkernel::Contacts instance
        // std::unique_ptr<meshkernel::CurvilinearGrid> m_curvilinearGrid; ///< Unique pointer to meshkernel::CurvilinearGrid instance

        // // Algorithms instances (interactivity)
        // std::unique_ptr<meshkernel::OrthogonalizationAndSmoothing> m_meshOrthogonalization;               ///< Unique pointer to meshkernel::OrthogonalizationAndSmoothing instance
        // std::unique_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines;             ///< Unique pointer to meshkernel::CurvilinearGridFromSplines instance
        // std::unique_ptr<meshkernel::CurvilinearGridOrthogonalization> m_curvilinearGridOrthogonalization; ///< Unique pointer to meshkernel::CurvilinearGridOrthogonalization instance
        // std::unique_ptr<meshkernel::CurvilinearGridLineShift> m_curvilinearGridLineShift;                 ///< Unique pointer to meshkernel::CurvilinearGridLineShift instance

        MeshKernelState(meshkernel::Projection projection) : m_projection(projection)
        {
            m_mesh1d = std::make_shared<meshkernel::Mesh1D>(projection);
            m_mesh2d = std::make_shared<meshkernel::Mesh2D>(projection);
            m_network1d = std::make_shared<meshkernel::Network1D>(projection);
            m_contacts = std::make_shared<meshkernel::Contacts>(*m_mesh1d, *m_mesh2d);
            m_curvilinearGrid = std::make_shared<meshkernel::CurvilinearGrid>(projection);
        }

        // Geometrical entities instances
        std::shared_ptr<meshkernel::Mesh1D> m_mesh1d;                   ///< Unique pointer to meshkernel::Mesh1D instance
        std::shared_ptr<meshkernel::Network1D> m_network1d;             ///< Unique pointer to meshkernel::Network1D instance
        std::shared_ptr<meshkernel::Mesh2D> m_mesh2d;                   ///< Unique pointer to meshkernel::Mesh2D instance
        std::shared_ptr<meshkernel::Contacts> m_contacts;               ///< Unique pointer to meshkernel::Contacts instance
        std::shared_ptr<meshkernel::CurvilinearGrid> m_curvilinearGrid; ///< Unique pointer to meshkernel::CurvilinearGrid instance

        // Algorithms instances (interactivity)
        std::shared_ptr<meshkernel::OrthogonalizationAndSmoothing> m_meshOrthogonalization;               ///< Unique pointer to meshkernel::OrthogonalizationAndSmoothing instance
        std::shared_ptr<meshkernel::CurvilinearGridFromSplines> m_curvilinearGridFromSplines;             ///< Unique pointer to meshkernel::CurvilinearGridFromSplines instance
        std::shared_ptr<meshkernel::CurvilinearGridOrthogonalization> m_curvilinearGridOrthogonalization; ///< Unique pointer to meshkernel::CurvilinearGridOrthogonalization instance
        std::shared_ptr<meshkernel::CurvilinearGridLineShift> m_curvilinearGridLineShift;                 ///< Unique pointer to meshkernel::CurvilinearGridLineShift instance

        // Exclusively owned state
        meshkernel::Projection m_projection{meshkernel::Projection::cartesian}; ///< Projection used by the meshes

        CurrentState m_state = CurrentState::Uninitialised;
    };


    static const std::string& toString(const MeshKernelState::CurrentState state);

} // namespace meshkernelapi

inline const std::string& meshkernelapi::toString(const MeshKernelState::CurrentState state)
{
    static std::string uninitialisedStr = "Uninitialised";
    static std::string validMeshStr = "ValidMesh";
    static std::string deletedMeshStr = "DeletedMesh";

    static std::string unknown = "UNKNOWN";

    using enum MeshKernelState::CurrentState;

    switch (state)
    {
    case Uninitialised:
        return uninitialisedStr;
    case ValidMesh:
        return validMeshStr;
    case DeletedMesh:
        return deletedMeshStr;
    default:
        return unknown;
    };
}
