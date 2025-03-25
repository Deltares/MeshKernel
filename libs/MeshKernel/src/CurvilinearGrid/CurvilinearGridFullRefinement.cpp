//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
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

#include "MeshKernel/CurvilinearGrid/CurvilinearGridFullRefinement.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridDeRefinement.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridRefinement.hpp"
#include "MeshKernel/CurvilinearGrid/UndoActions/CurvilinearGridRefinementUndoAction.hpp"
#include "MeshKernel/Exceptions.hpp"

meshkernel::UndoActionPtr meshkernel::CurvilinearGridFullRefinement::Compute(CurvilinearGrid& grid,
                                                                             const int mRefinement,
                                                                             const int nRefinement) const
{
    if (mRefinement == 0 && nRefinement == 0)
    {
        throw ConstraintError("Incorrect refinement value. One of the two values must be larger than 0",
                              mRefinement, nRefinement);
    }

    if (mRefinement == constants::missing::intValue || nRefinement == constants::missing::intValue)
    {
        throw ConstraintError("Incorrect refinement value. Value should not be the invalid value: m-refinement = {}, n-refinement = {}",
                              mRefinement, nRefinement);
    }

    if (mRefinement == 1 && nRefinement == 1)
    {
        // nothing to do
        return nullptr;
    }

    if (mRefinement > 1 && nRefinement > 1)
    {
        return Compute(grid, mRefinement, nRefinement);
    }

    // Estimate the dimension of the refined grid
    const UInt maxM = grid.NumM();
    const auto left = CurvilinearGridNodeIndices(0, 0);
    const auto right = CurvilinearGridNodeIndices(0, maxM);

    std::unique_ptr<CompoundUndoAction> undoAction = CompoundUndoAction::Create();
    if (mRefinement > 1)
    {
        CurvilinearGridRefinement curvilinearGridRefinement(grid, mRefinement);
        curvilinearGridRefinement.SetBlock(left, right);
        undoAction->Add(curvilinearGridRefinement.Compute());
    }
    else if (mRefinement < 1)
    {
        CurvilinearGridDeRefinement curvilinearGridDeRefinement(grid, -mRefinement);
        curvilinearGridDeRefinement.SetBlock(left, right);
        undoAction->Add(curvilinearGridDeRefinement.Compute());
    }

    const UInt maxN = grid.NumN();
    const auto bottom = CurvilinearGridNodeIndices(0, 0);
    const auto top = CurvilinearGridNodeIndices(maxN, 0);
    if (nRefinement > 1)
    {
        CurvilinearGridRefinement curvilinearGridRefinement(grid, nRefinement);
        curvilinearGridRefinement.SetBlock(bottom, top);
        undoAction->Add(curvilinearGridRefinement.Compute());
    }
    else if (nRefinement < 1)
    {
        CurvilinearGridDeRefinement curvilinearGridDeRefinement(grid, -nRefinement);
        curvilinearGridDeRefinement.SetBlock(bottom, top);
        undoAction->Add(curvilinearGridDeRefinement.Compute());
    }

    return undoAction;
}

meshkernel::UndoActionPtr meshkernel::CurvilinearGridFullRefinement::Compute(CurvilinearGrid& grid,
                                                                             const UInt mRefinement,
                                                                             const UInt nRefinement) const
{

    if (mRefinement == 0 || nRefinement == 0)
    {
        throw ConstraintError("Incorrect refinement value. Value must be greater than 0: m-refinement = {}, n-refinement = {}",
                              mRefinement, nRefinement);
    }

    if (mRefinement == constants::missing::uintValue || nRefinement == constants::missing::uintValue)
    {
        throw ConstraintError("Incorrect refinement value. Value should not be the invalid value: m-refinement = {}, n-refinement = {}",
                              mRefinement, nRefinement);
    }

    if (mRefinement == 1 && nRefinement == 1)
    {
        // nothing to do
        return nullptr;
    }

    Splines splines(grid);

    // Estimate the dimension of the refined grid
    const UInt maxM = (grid.NumM() - 1) * mRefinement + 1;
    const UInt maxN = (grid.NumN() - 1) * nRefinement + 1;

    // Local vector for each curvilinear grid face
    std::vector<Point> bottomRefinement(mRefinement + 1);
    std::vector<Point> topRefinement(mRefinement + 1);
    std::vector<Point> leftRefinement(nRefinement + 1);
    std::vector<Point> rightRefinement(nRefinement + 1);

    // The nodes of the refined grid
    lin_alg::Matrix<Point> refinedGrid(maxN, maxM);

    UInt refinedN = 0;

    for (UInt currentN = 0; currentN < grid.NumN() - 1; ++currentN)
    {
        UInt refinedM = 0;

        for (UInt currentM = 0; currentM < grid.NumM() - 1; ++currentM)
        {
            // Only if all grid nodes of the face are valid, perform transfinite interpolation
            if (ValidFace(grid, currentM, currentN))
            {
                ComputeRefinedElementEdges(splines,
                                           grid.NumM() + currentN,
                                           currentM,
                                           mRefinement,
                                           bottomRefinement,
                                           topRefinement);

                ComputeRefinedElementEdges(splines,
                                           currentM,
                                           currentN,
                                           nRefinement,
                                           leftRefinement,
                                           rightRefinement);

                // Perform transfinite interpolation on the current curvilinear face
                const auto localGrid = DiscretizeTransfinite(bottomRefinement,
                                                             topRefinement,
                                                             leftRefinement,
                                                             rightRefinement,
                                                             grid.projection(),
                                                             mRefinement,
                                                             nRefinement);

                // Copy the local grid into the refined grid block
                refinedGrid.block(refinedN, refinedM, nRefinement + 1, mRefinement + 1) = localGrid;
            }

            refinedM += mRefinement;
        }

        refinedN += nRefinement;
    }

    std::unique_ptr<CurvilinearGridRefinementUndoAction> undoAction = CurvilinearGridRefinementUndoAction::Create(grid);

    // Substitute original grid with the refined one
    grid.SetGridNodes(refinedGrid);

    return undoAction;
}

bool meshkernel::CurvilinearGridFullRefinement::ValidFace(const CurvilinearGrid& grid,
                                                          const UInt m,
                                                          const UInt n) const
{
    return grid.GetNode(n, m).IsValid() &&
           grid.GetNode(n + 1, m).IsValid() &&
           grid.GetNode(n, m + 1).IsValid() &&
           grid.GetNode(n + 1, m + 1).IsValid();
}

void meshkernel::CurvilinearGridFullRefinement::ComputeRefinedElementEdges(const Splines& splines,
                                                                           const UInt splineIndex,
                                                                           const UInt currentIndex,
                                                                           const UInt refinement,
                                                                           std::vector<Point>& elementSide1,
                                                                           std::vector<Point>& elementSide2) const
{
    elementSide1.resize(refinement + 1);
    elementSide2.resize(refinement + 1);

    for (UInt i = 0; i < refinement + 1; ++i)
    {
        const auto interpolationPoint = static_cast<double>(currentIndex) + static_cast<double>(i) / static_cast<double>(refinement);
        elementSide1[i] = splines.Evaluate(splineIndex, interpolationPoint);
        elementSide2[i] = splines.Evaluate(splineIndex + 1, interpolationPoint);
    }
}
