#include "MeshKernel/CurvilinearGrid/CurvilinearGridEntireRefinement.hpp"
#include "MeshKernel/Splines.hpp"

void meshkernel::CurvilinearGridEntireRefinement::Compute(CurvilinearGrid& grid, UInt mRefinement, UInt nRefinement) const
{
    Splines splines (grid);

    // Estimate the dimension of the refined grid
    const UInt maxM = (grid.NumM() - 1) * mRefinement + 1;
    const UInt maxN = (grid.NumN() - 1) * nRefinement + 1;

    std::cout << "new grid size (n,m): " << maxN << "  "<< maxM << " refinement levels (m,n): " << mRefinement << "  " << nRefinement << std::endl;

    // Local vector for each curvilinear grid face
    std::vector<Point> bottomRefinement(mRefinement);
    std::vector<Point> topRefinement(mRefinement);
    std::vector<Point> leftRefinement(nRefinement);
    std::vector<Point> rightRefinement(nRefinement);

    // The refined grid
    lin_alg::Matrix<Point> refinedGrid(maxN, maxM);

    UInt refinedN = 0;

    for (UInt currentN = 0; currentN < grid.NumN() - 1; ++currentN)
    {
        UInt refinedM = 0;

        for (UInt currentM = 0; currentM < grid.NumM() - 1; ++currentM)
        {

            // Only if all grid nodes of the face are valid, perform transfinite interpolation
            if (grid.GetNode(currentN, currentM).IsValid() &&
                grid.GetNode(currentN + 1, currentM).IsValid() &&
                grid.GetNode(currentN, currentM + 1).IsValid() &&
                grid.GetNode(currentN + 1, currentM + 1).IsValid())
            {
                // Calculate m-direction spline points
                bottomRefinement.clear();
                topRefinement.clear();

                // Calculate n-direction spline points
                leftRefinement.clear();
                rightRefinement.clear();

                // for (UInt n = 0; n < mRefinement + 1; ++n)
                // {
                //     const auto splineIndex = grid.NumM() + currentN;
                //     const auto interpolationPoint = static_cast<double>(currentM) + static_cast<double>(n) / static_cast<double>(mRefinement);
                //     bottomRefinement.emplace_back(splines.Evaluate(splineIndex, interpolationPoint));
                //     topRefinement.emplace_back(splines.Evaluate(splineIndex + 1, interpolationPoint));
                // }

                // for (UInt m = 0; m < nRefinement + 1; ++m)
                // {
                //     const auto splineIndex = currentM;
                //     const auto interpolationPoint = static_cast<double>(currentN) + static_cast<double>(m) / static_cast<double>(nRefinement);
                //     leftRefinement.emplace_back(splines.Evaluate(splineIndex, interpolationPoint));
                //     rightRefinement.emplace_back(splines.Evaluate(splineIndex + 1, interpolationPoint));
                // }

                for (UInt n = 0; n < nRefinement + 1; ++n)
                {
                    const auto splineIndex = currentM;
                    const auto interpolationPoint = static_cast<double>(currentN) + static_cast<double>(n) / static_cast<double>(nRefinement);
                    bottomRefinement.emplace_back(splines.Evaluate(splineIndex, interpolationPoint));
                    topRefinement.emplace_back(splines.Evaluate(splineIndex + 1, interpolationPoint));
                }

                for (UInt m = 0; m < mRefinement + 1; ++m)
                {
                    const auto splineIndex = grid.NumM() + currentN;
                    const auto interpolationPoint = static_cast<double>(currentM) + static_cast<double>(m) / static_cast<double>(mRefinement);
                    leftRefinement.emplace_back(splines.Evaluate(splineIndex, interpolationPoint));
                    rightRefinement.emplace_back(splines.Evaluate(splineIndex + 1, interpolationPoint));
                }


                if (currentM == 0)
                {
                    std::cout << "bottom ref: ";

                    for (UInt i = 0; i < bottomRefinement.size (); ++i)
                    {
                        std::cout << "{" << bottomRefinement[i].x << ", " << bottomRefinement[i].y << " } -- ";
                    }

                    std::cout << std::endl;
                    std::cout << "left ref: ";

                    for (UInt i = 0; i < leftRefinement.size (); ++i)
                    {
                        std::cout << "{" << leftRefinement[i].x << ", " << leftRefinement[i].y << " } -- ";
                    }

                    std::cout << std::endl;
                }

                // Perform transfinite interpolation on the current curvilinear face
                const auto localGrid = DiscretizeTransfinite(bottomRefinement,
                                                             topRefinement,
                                                             leftRefinement,
                                                             rightRefinement,
                                                             grid.projection(),
                                                             mRefinement,
                                                             nRefinement);


                std::cout << " -------------------------------- " << std::endl
                          << " element refinement " << localGrid.rows () << "  " << localGrid.cols () << std::endl;

                // Copy the local grid into the refined grid
                for (Eigen::Index n = 0; n < nRefinement + 1; ++n)
                {
                    for (Eigen::Index m = 0; m < mRefinement + 1; ++m)
                    {
                        std::cout << localGrid(n, m).x << ", " << localGrid(n, m).y << std::endl;
                        refinedGrid(refinedN + n, refinedM + m) = localGrid(n, m);
                    }
                }

                std::cout << std::endl;

            }

            refinedM += mRefinement;
        }

        refinedN += nRefinement;
    }

    std::cout << "old grid size: " << grid.NumN () << "  "<< grid.NumM () << std::endl;
    // Substitute original grid with the refined one
    grid.SetGridNodes(refinedGrid);

    std::cout << "new grid size: " << grid.NumN () << "  "<< grid.NumM () << std::endl;

}
