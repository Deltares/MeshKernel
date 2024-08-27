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

#include "MeshKernelApi/MKStateUndoAction.hpp"

std::unique_ptr<meshkernelapi::MKStateUndoAction> meshkernelapi::MKStateUndoAction::Create(MeshKernelState& mkState)
{
    return std::make_unique<MKStateUndoAction>(mkState);
}

meshkernelapi::MKStateUndoAction::MKStateUndoAction(MeshKernelState& mkState) : m_mkStateReference(mkState)
{
    // Copy the pointers
    m_mkState.m_mesh1d = mkState.m_mesh1d;
    m_mkState.m_network1d = mkState.m_network1d;
    m_mkState.m_mesh2d = mkState.m_mesh2d;
    m_mkState.m_contacts = mkState.m_contacts;
    m_mkState.m_curvilinearGrid = mkState.m_curvilinearGrid;
    m_mkState.m_meshOrthogonalization = mkState.m_meshOrthogonalization;
    m_mkState.m_curvilinearGridFromSplines = mkState.m_curvilinearGridFromSplines;
    m_mkState.m_curvilinearGridOrthogonalization = mkState.m_curvilinearGridOrthogonalization;
    m_mkState.m_curvilinearGridLineShift = mkState.m_curvilinearGridLineShift;
    m_mkState.m_projection = mkState.m_projection;
}

void meshkernelapi::MKStateUndoAction::SwapContents()
{
    std::swap(m_mkState.m_mesh1d, m_mkStateReference.m_mesh1d);
    std::swap(m_mkState.m_mesh2d, m_mkStateReference.m_mesh2d);
    std::swap(m_mkState.m_network1d, m_mkStateReference.m_network1d);
    std::swap(m_mkState.m_contacts, m_mkStateReference.m_contacts);
    std::swap(m_mkState.m_curvilinearGrid, m_mkStateReference.m_curvilinearGrid);
}

void meshkernelapi::MKStateUndoAction::DoCommit()
{
    SwapContents();
}

void meshkernelapi::MKStateUndoAction::DoRestore()
{
    SwapContents();
}
