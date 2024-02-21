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

#include "MeshKernel/UndoAction.hpp"

namespace meshkernel
{

    template <typename DerivedUndoAction, class Mesh>
    class BaseMeshUndoAction : public UndoAction
    {
    public:
        BaseMeshUndoAction(Mesh& mesh) : mesh_(mesh) {}

    private:
        void DoCommit();

        void DoRestore();

        Mesh& mesh_;
    };

    // template <typename DerivedUndoAction, class Mesh>
    // class BaseMeshUndoAction : public UndoAction
    // {
    // public:
    //     BaseMeshUndoAction(Mesh& mesh) : mesh_(mesh) {}

    // private:
    //     void doCommit();

    //     void doRestore();

    //     Mesh& mesh_;
    // };

} // namespace meshkernel

template <typename DerivedUndoAction, class Mesh>
void meshkernel::BaseMeshUndoAction<DerivedUndoAction, Mesh>::DoCommit()
{
    mesh_.commit(*static_cast<DerivedUndoAction*>(this));
    // mesh_.commit(static_cast<DerivedUndoAction&>(*this));
}

template <typename DerivedUndoAction, class Mesh>
void meshkernel::BaseMeshUndoAction<DerivedUndoAction, Mesh>::DoRestore()
{
    mesh_.restore(*static_cast<DerivedUndoAction*>(this));
    // mesh_.restore(static_cast<DerivedUndoAction&>(*this));
}
