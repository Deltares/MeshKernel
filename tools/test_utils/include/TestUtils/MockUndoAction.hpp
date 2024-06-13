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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MeshKernel/UndoActions/UndoAction.hpp"

/// @brief Mock class for UndoAction's
class MockUndoAction : public meshkernel::UndoAction
{
public:
    /// @brief Mock the DoCommit function
    MOCK_METHOD(void, DoCommit, (), (override));

    /// @brief Mock the DoRestore function
    MOCK_METHOD(void, DoRestore, (), (override));
};
