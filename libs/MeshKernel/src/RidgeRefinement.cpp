#include "MeshKernel/RidgeRefinement.hpp"

meshkernel::Hessian::Hessian(const UInt dim1, const UInt dim2, const UInt dim3)
{
    resize(dim1, dim2, dim3);
}

void meshkernel::Hessian::resize(const UInt dim1, const UInt dim2, const UInt dim3)
{
    m_hessian.resize(dim1);

    for (UInt i = 0; i < dim1; ++i)
    {
        m_hessian[i].resize(dim2, dim3);
    }

    zero();
}

void meshkernel::Hessian::zero()
{
    for (UInt i = 0; i < m_hessian.size(); ++i)
    {
        m_hessian[i].setZero();
    }
}

void meshkernel::RidgeRefinement::GetValidPoints(const std::vector<Point>& samplePoints, std::vector<Point>& validSamplePoints, std::vector<UInt>& iperm) const
{
    validSamplePoints.resize(samplePoints.size());
    size_t validCount = 0;

    // Do not copy invalid points.
    for (size_t i = 0; i < samplePoints.size(); ++i)
    {
        if (samplePoints[i].IsValid())
        {
            validSamplePoints[validCount] = samplePoints[i];
            iperm[validCount] = i;
            ++validCount;
        }
    }

    validSamplePoints.resize(validCount);
}

void meshkernel::RidgeRefinement::RemoveDuplicates(std::vector<Point>& samplePoints, std::vector<double>& sampleData, std::vector<UInt>& sampleIndices /*ipsam*/) const
{
    std::vector<Point> validSamplePoints;
    // change name and formal parameter name
    std::vector<UInt> iperm(samplePoints.size(), constants::missing::uintValue);
    std::vector<UInt> newNode(samplePoints.size(), constants::missing::uintValue);
    GetValidPoints(samplePoints, validSamplePoints, iperm);

    RTree kdTree;

    kdTree.BuildTree(validSamplePoints);

    std::cout << " kdtree size  " << kdTree.Size() << std::endl;

    bool duplicatesFound = true;
    UInt nummerged = 0;
    UInt sampleCount = samplePoints.size();

    double tolerance = 1.0e-4;
    double searchRadius = tolerance * tolerance;

    std::cout.precision(10);

    while (duplicatesFound)
    {
        std::cout << "Duplicates found" << std::endl;
        duplicatesFound = false;

        // Check for duplicate sample data points
        for (size_t ii = 0; ii < sampleCount; ++ii)
        {
            if (iperm[ii] == constants::missing::uintValue)
            {
                // already merged
                continue;
            }

            UInt i = iperm[ii];
            kdTree.SearchPoints(samplePoints[iperm[i]], searchRadius);

            if (kdTree.GetQueryResultSize() > 1)
            {
                for (UInt k = 0; k < kdTree.GetQueryResultSize(); ++k)
                {
                    UInt jj = kdTree.GetQueryResult(k);
                    UInt j = iperm[jj];

                    // Which point comparison to use == or IsEqual (wth a tolerance)
                    if (j != i && j != constants::missing::uintValue && IsEqual(samplePoints[j], samplePoints[i], tolerance))
                    {
                        iperm[jj] = constants::missing::uintValue;
                        samplePoints[j].SetInvalid();
                        duplicatesFound = true;
                        ++nummerged;
                    }
                }
            }
        }

        UInt count = 0;

        // Removed duplicates
        for (UInt i = 0; i < sampleCount; ++i)
        {
            if (samplePoints[i].IsValid() && sampleData[i] != constants::missing::doubleValue)
            {
                samplePoints[count] = samplePoints[i];
                sampleData[count] = sampleData[i];
                newNode[i] = count;
                ++count;
            }
        }

        UInt k = 0;

        // Update permutation
        for (UInt i = 0; i < sampleCount; ++i)
        {
            UInt j = sampleIndices[i];
            if (newNode[j] != constants::missing::uintValue)
            {
                sampleIndices[k] = newNode[j];
                ++k;
            }
        }

        sampleCount = k;
    }

    UInt newSize = samplePoints.size() - nummerged;

    samplePoints.resize(newSize);
    sampleData.resize(newSize);
    sampleIndices.resize(newSize);
}

void meshkernel::RidgeRefinement::TidySamples(std::vector<Point>& samplePoints, std::vector<double>& sampleData) const
{
    std::vector<UInt> sampleIndices(samplePoints.size());
    std::iota(sampleIndices.begin(), sampleIndices.end(), 0);

    auto compareX = [samplePoints](const size_t i, const size_t j)
    { return samplePoints[i].x < samplePoints[j].x; };

    std::ranges::sort(sampleIndices, compareX);

    RemoveDuplicates(samplePoints, sampleData, sampleIndices);
}

void meshkernel::RidgeRefinement::smoothSamples(const std::vector<Point>& samplePoints,
                                                const std::vector<double>& sampleData,
                                                Hessian& hessian) const
{

    lin_alg::MatrixColMajor<double> zsdum(mesh.m_numM, mesh.m_numN);
    zsdum.setZero();

    // May be better to reorder the loop (j, then i), can then use a count for the sampleData.
    for (UInt i = 0; j < hessian.size(1); ++i)
    {
        for (UInt j = 0; j < hessian.size(2); ++j)
        {
            hessian(0, i, j) = sampleData[i + hessian.size(1) * j];
        }
    }

    for (int iter = 1; iter <= numberOfSmoothingIterations; ++iter)
    {
        double af = static_cast<double>(iter - 1) / static_cast<double>(std::max(numberOfSmoothingIterations - 1, 1));
        zsdum = hessian.getMatrix(0);

        for (UInt j = 1; j < hessian.size(2) - 1; ++j)
        {
            for (UInt i = 1; j < hessian.size(1) - 1; ++i)
            {

                if (zsdum(i, j) == constants::missing::doubleValue)
                {
                    continue;
                }

                // Needs tidying up, seems not to be really necessary
                ciL = 1.0;
                ciR = 1.0;
                cjL = 1.0;
                cjR = 1.0;

                if (zsdum(i - 1, j) == constants::missing::doubleValue)
                {
                    ciL = 0.0;
                }
                if (zsdum(i + 1, j) == constants::missing::doubleValue)
                {
                    ciR = 0.0;
                }
                if (zsdum(i, j - 1) == constants::missing::doubleValue)
                {
                    cjL = 0.0;
                }
                if (zsdum(i, j + 1) == constants::missing::doubleValue)
                {
                    cjR = 0.0;
                }

                if (ciL * ciR * cjL * cjR == 0.0)
                {
                    continue;
                }

                c0 = ciL + ciR + cjL + cjR;

                if (std::abs(c0) < 0.5)
                {
                    continue;
                }

                zss(1, i, j) = (1.0 - sigma) * zsdum(i, j) +
                               sigma * (ciL * zsdum(i - 1, j) + ciR * zsdum(i + 1, j) + cjL * zsdum(i, j - 1) + ckR * zsdum(i, j + 1)) / c0;
            }
        }
    }
}

void meshkernel::RidgeRefinement::computeGradient(const std::vector<Point>& samplePoints,
                                                  const std::vector<double>& sampleData,
                                                  const Hessian& hessian,
                                                  const UInt ip0,
                                                  const UInt ip1,
                                                  const UInt ip0L,
                                                  const UInt ip0R,
                                                  const UInt ip1L,
                                                  const UInt ip1R,
                                                  Eigen::Vector2d& gradient,
                                                  Eigen::Vector2d& S) const
{
    //   compute the gradient in a control volume defined by the polygon (0-R-1-L)
    //
    //       0L-----1L
    //        |     |
    //        |  L  |
    //        | / \ |
    //        |/   \|
    //        0-----1
    //        |\   /|
    //        | \ / |
    //        |  R  |
    //        |     |
    //       0R-----1R
    //
    //  0 and 1 are sample points
    //  L and R are interpolated at sample cell centers

    gradient(0) = DMISS;
    gradient(1) = DMISS;
    S(0) = 0.0;
    S(1) = 0.0;
    dareaL = 0.0;
    dareaR = 0.0;

    x0 = xs(ip0);
    y0 = ys(ip0);
    z0 = zss(1, ip0);

    x1 = xs(ip1);
    y1 = ys(ip1);
    z1 = zss(1, ip1);

    if (x0.eq.DMISS.or.y1.eq.DMISS.or.x1.eq.DMISS.or.y1.eq.DMISS)
    {
        return;
    }

    leftPoint = 0.25 * (samplePoints[ip0] + samplePoints[ip1] + samplePoints[ip0L] + samplePoints[ip1L]);
    rightPoint = 0.25 * (samplePoints[ip0] + samplePoints[ip1] + samplePoints[ip0R] + samplePoints[ip1R]);

    cx1 = -GetDelta(leftPoint, rightPoint, projection);
    cxL = -GetDelta(samplePoints[ip0], samplePoints[ip1], projection);

    cx0 = -cx1;
    cxR = -cxL;

    darea = 0.5 * (dot(cx0, x0) + dot(cx1, x1) + dot(cxL, xL) + dot(cxR, xR));

    // darea = 0.5 * (cx0.x () * x0 + cy0.);
    // darea = 0.5 * (cx0 * x0 + cy0 * y0 + cx1 * x1 + cy1 * y1 + cxL * xL + cyL * yL + cxR * xR + cyR * yR);

    // !     gradx and grady can be composed

    if (zss(1, ip0) != DMISS && zss(1, ip1) != DMISS && &zss(1, ip0L) != DMISS && zss(1, ip0R) != DMISS && &zss(1, ip1L) != DMISS && zss(1, ip1R) != DMISS)
    {
        double zL = 0.25 * (zss(1, ip0) + zss(1, ip1) + zss(1, ip0L) + zss(1, ip1L));
        double zR = 0.25 * (zss(1, ip0) + zss(1, ip1) + zss(1, ip0R) + zss(1, ip1R));
        gradient(0) = (cx1.x() * z1 + cxL.x() * zL + cx0.x() * z0 + cxR.x() * zR) / darea;
        gradient(1) = (cy1.y() * z1 + cyL.y() * zL + cy0.y() * z0 + cyR.y() * zR) / darea;
    }

    S = 2.0 * cx1;
    drealL = 0.5 * std::abs(OuterProductTwoSegments(x0, xR, x0, xL, projection));
    drealR = 0.5 * std::abs(OuterProductTwoSegments(x1, xR, x1, xL, projection));
}

void meshkernel::RidgeRefinement::computeSampleGradient(const std::vector<Point>& samplePoints,
                                                        const std::vector<double>& sampleData,
                                                        const UInt direction,
                                                        const UInt i,
                                                        const UInt j,
                                                        Eigen::Vector2d& gradient,
                                                        Eigen::Vector2d& sn,
                                                        double& dareaL,
                                                        double& dereaR) const
{
    gradient(0) = 0.0;
    gradient(1) = 0.0;
    sn(0) = 0.0;
    sn(1) = 0.0;

    dareaL = 0.0;
    dareaR = 0.0;

    if (direction == 0)
    {
        //      i-edge gradient at (i+1/2,j) location
        //         control volume:
        //
        //                   L:(i+1/2,j+1/2)
        //                  / \
        //                 /   \
        //          0:(i,j)-----1:(i+1,j)
        //                 \   /
        //                  \ /
        //                   R:(i+1/2,j-1/2)

        ip0 = i + hessian.size(1) * (j - 1);                           // ! pointer to (i,j)
        ip1 = i + 1 + hessian.size(1) * (j - 1);                       // ! pointer to (i+1,j)
        ip0L = i + hessian.size(1) * (std::min(j + 1, MYSAM) - 1);     // ! pointer to (i,j+1)
        ip0R = i + hessian.size(1) * (std::max(j - 1, 1) - 1);         // ! pointer to (i,j-1)
        ip1L = i + 1 + hessian.size(1) * (std::min(j + 1, MYSAM) - 1); // ! pointer to (i+1,j+1)
        ip1R = i + 1 + hessian.size(1) * (std::max(j - 1, 1) - 1);     // ! pointer to (i+1,j-1)
        computeGradient(samplePoints, sampleData, hessian, ip0, ip1, ip0L, ip0R, ip1L, ip1R, gradient, S, DareaL, DareaR);
    }
    else if (direction == 1)
    {
        //      j-edge gradient at (i,j+1/2) location
        //        control volume:
        //
        //                   1:(i,j+1)
        //                  / \
        //                 /   \
        //  L:(i-1/2,j+1/2)-----R:(i+1/2,j+1/2)
        //                 \   /
        //                  \ /
        //                   0:(i,j)

        ip0 = i + hessian.size(1) * (j - 1);                                 //              ! pointer to (i,j)
        ip1 = i + hessian.size(1) * (j);                                     //              ! pointer to (i,j+1)
        ip0L = std::max(i - 1, 1) + hessian.size(1) * (j - 1);               //              ! pointer to (i-1,j)
        ip0R = std::min(i + 1, hessian.size(1)) + hessian.size(1) * (j - 1); //              ! pointer to (i+1,j)
        ip1L = std::max(i - 1, 1) + hessian.size(1) * (j);                   //              ! pointer to (i-1,j+1)
        ip1R = std::min(i + 1, hessian.size(1)) + hessian.size(1) * (j);     //              ! pointer to (i+1,j+1)
        computeGradient(samplePoints, sampleData, hessian, ip0, ip1, ip0L, ip0R, ip1L, ip1R, gradient, S, DareaL, DareaR);
    }
}

void meshkernel::RidgeRefinement::computeHessian(const std::vector<Point>& samplePoints,
                                                 const std::vector<double>& sampleData,
                                                 Hessian& hessian) const
{
    double dh = std::min(ComputeDistance(samplePoints[0], samplePoints[1], projection),
                         ComputeDistance(samplePoints[hessian.size(1) - 1], samplePoints[hessian.size(1) - 1], projection));

    if (hessian.size(1) < 3 || hessian.size(2) < 2)
    {
        return;
    }

    for (UInt i = 1; i < hessian.size(1) - 1; ++i)
    {
        double af = static_cast<double>(i - 1) / static_cast<double>(std::max(hessian.size(1) - 3, 1));

        for (UInt j = 1; j < hessian.size(2) - 1; ++j)
        {
            ip = i + (j - 1) * hessian.size(1);

            if ((std::abs(samplePoints[ip].x - 87270.0) < 1.0e-8) && (std::abs(samplePoints[ip].y - 415570.0) < 1.0e-8))
            {
                continue;
            }

            double zxx = 0.0;
            double zxy = 0.0;
            double zyx = 0.0;
            double zyy = 0.0;

            double UU = 0.0;
            double VV = 0.0;

            double zx = 0.0;
            double zy = 0.0;

            double S = 0.0;
            UInt k = 0;
            bool hasRidge = false;

            for ()
            {
                computeSampleGradient(samplePoints, sampleData, 0, i, j, gradientiR, sniR, dareaiR, dum);

                if (gradientiR(0) == constants::missing::doubleValue)
                {
                    break;
                }

                computeSampleGradient(samplePoints, sampleData, 1, i - 1, j, gradientiL, sniL, dum, dareaiL);

                if (gradientiL(0) == constants::missing::doubleValue)
                {
                    break;
                }

                computeSampleGradient(samplePoints, sampleData, 1, i, j, gradientjR, snij, dareajR, dum);

                if (gradientjR(0) == constants::missing::doubleValue)
                {
                    break;
                }

                computeSampleGradient(samplePoints, sampleData, 1, i, j - 1, gradientjL, snjL, dum, dareajL);

                if (gradientjL(0) == constants::missing::doubleValue)
                {
                    break;
                }

                area = dareaiL + dareaiR + dareajL + dareajR;

                zx = (0.5d0 * (zss(1, i + 1, j) + zss(1, i, j)) * SniR(1) - 0.5d0 * (zss(1, i - 1, j) + zss(1, i, j)) * SniL(1) +
                      0.5d0 * (zss(1, i, j + 1) + zss(1, i, j)) * SnjR(1) - 0.5d0 * (zss(1, i, j - 1) + zss(1, i, j)) * SnjL(1)) /
                     area;
                zy = (0.5d0 * (zss(1, i + 1, j) + zss(1, i, j)) * SniR(2) - 0.5d0 * (zss(1, i - 1, j) + zss(1, i, j)) * SniL(2) +
                      0.5d0 * (zss(1, i, j + 1) + zss(1, i, j)) * SnjR(2) - 0.5d0 * (zss(1, i, j - 1) + zss(1, i, j)) * SnjL(2)) /
                     area;

                VV(0, 0) = (gradiR(1) * SniR(1) - gradiL(1) * SniL(1) + gradjR(1) * SnjR(1) - gradjL(1) * SnjL(1)) / area;
                VV(0, 1) = (gradiR(1) * SniR(2) - gradiL(1) * SniL(2) + gradjR(1) * SnjR(2) - gradjL(1) * SnjL(2)) / area;
                VV(1, 0) = (gradiR(2) * SniR(1) - gradiL(2) * SniL(1) + gradjR(2) * SnjR(1) - gradjL(2) * SnjL(1)) / area;
                VV(1, 1) = (gradiR(2) * SniR(2) - gradiL(2) * SniL(2) + gradjR(2) * SnjR(2) - gradjL(2) * SnjL(2)) / area;

                // Eigendecompostion
                Eigen::EigenSolver<Matrix2d> eigensolver(VV);
                jacobi(VV, 2, 2, S, UU, nrot);

                k = std::abs(S(1)) > std::abs(S(2)) ? 0 : 1;

                hessian(2, i, j) = VV(0, k);
                hessian(3, i, j) = VV(1, k);
                hessian(4, i, j) = S(k) * area;
                hessian(5, i, j) = -(VV(1, k) * zx + VV(2, k) * zy) / (S(k) + 1.0e-8);
            }
        }
    }
}

void meshkernel::RidgeRefinement::prepareSampleForHessian(const std::vector<Point>& samplePoints,
                                                          const std::vector<double>& sampleData,
                                                          Hessian& hessian) const
{
    smoothSamples(sampleData, numberOfSmoothingIterations, hessian);
    computeHessian(sampleData, hessian);
}

void meshkernel::RidgeRefinement::Compute(Mesh2D& mesh, const std::vector<Point>& rawSamplePoints, const std::vector<double>& rawSampleData) const
{
    (void)mesh;
    std::vector<Point> samplePoints(rawSamplePoints);
    std::vector<double> sampleData(rawSampleData);

    TidySamples(samplePoints, sampleData);
    BoundingBox boundingBox(samplePoints);

    RTree kdTree;
    kdTree.BuildTree(samplePoints);

    // Hessian hessian(5, mesh.);
    // prepareSampleForHessian(samplePoints, sampleData);
}

void meshkernel::RidgeRefinement::Compute(const CurvilinearMesh& mesh, const std::vector<Point>& rawSamplePoints, const std::vector<double>& rawSampleData) const
{
    (void)mesh;
    std::vector<Point> samplePoints(rawSamplePoints);
    std::vector<double> sampleData(rawSampleData);

    TidySamples(samplePoints, sampleData);
    BoundingBox boundingBox(samplePoints);

    RTree kdTree;
    kdTree.BuildTree(samplePoints);

    Hessian hessian(5, mesh.m_numM, mesh.m_numN);
    prepareSampleForHessian(mesh, samplePoints, sampleData, hessian);
}
