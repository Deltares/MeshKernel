#include "MeshKernel/HessianCalculator.hpp"
#include "MeshKernel/Operations.hpp"

#include <Eigen/Core>

void meshkernel::HessianCalculator::SmoothSamples(const std::vector<Sample>& sampleData,
                                                  const UInt numberOfSmoothingIterations,
                                                  Hessian& hessian)
{

    const double sigma = 0.5;

    // Check the dimension are correct.
    MatrixColMajor zsdum(hessian.size(1), hessian.size(2));
    zsdum.setZero();

    [[maybe_unused]] UInt count = 0;

    // May be better to reorder the loop (j, then i), can then use a count for the sampleData.
    for (UInt i = 0; i < hessian.size(1); ++i)
    {
        for (UInt j = 0; j < hessian.size(2); ++j)
        {
            // hessian(0, i, j) = sampleData[count].value;
            hessian(0, i, j) = sampleData[i + hessian.size(1) * j].value;
            ++count;
        }
    }

    for (UInt iter = 1; iter <= numberOfSmoothingIterations; ++iter)
    {
        // double af = static_cast<double>(iter - 1) / static_cast<double>(std::max(numberOfSmoothingIterations - 1, 1));
        zsdum = hessian.getMatrix(0);

        for (UInt j = 1; j < hessian.size(2) - 1; ++j)
        {
            for (UInt i = 1; i < hessian.size(1) - 1; ++i)
            {

                if (zsdum(i, j) == constants::missing::doubleValue)
                {
                    continue;
                }

                // Needs tidying up, seems not to be really necessary
                double ciL = 1.0;
                double ciR = 1.0;
                double cjL = 1.0;
                double cjR = 1.0;

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

                double c0 = ciL + ciR + cjL + cjR;

                if (std::abs(c0) < 0.5)
                {
                    continue;
                }

                hessian(0, i, j) = (1.0 - sigma) * zsdum(i, j) +
                                   sigma * (ciL * zsdum(i - 1, j) + ciR * zsdum(i + 1, j) + cjL * zsdum(i, j - 1) + cjR * zsdum(i, j + 1)) / c0;
            }
        }
    }
}

void meshkernel::HessianCalculator::ComputeGradient(const std::vector<Sample>& samplePoints,
                                                    const Projection projection,
                                                    const Hessian& hessian,
                                                    const UInt ip0,
                                                    const UInt ip1,
                                                    const UInt ip0L,
                                                    const UInt ip0R,
                                                    const UInt ip1L,
                                                    const UInt ip1R,
                                                    meshkernel::Vector& gradient,
                                                    meshkernel::Vector& S,
                                                    double& dareaL,
                                                    double& dareaR)
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

    gradient[0] = constants::missing::doubleValue;
    gradient[1] = constants::missing::doubleValue;
    S[0] = 0.0;
    S[1] = 0.0;
    dareaL = 0.0;
    dareaR = 0.0;

    Point x0 = samplePoints[ip0];
    double z0 = samplePoints[ip0].value;

    const MatrixColMajor& zss = hessian.getMatrix(0);

    Point x1 = samplePoints[ip1];
    double z1 = zss(ip1); // zss(1, ip1);

    if (x0.x == constants::missing::doubleValue ||
        x0.y == constants::missing::doubleValue ||
        x1.x == constants::missing::doubleValue ||
        x1.y == constants::missing::doubleValue)
    {
        return;
    }

    Point leftPoint = 0.25 * (samplePoints[ip0] + samplePoints[ip1] + samplePoints[ip0L] + samplePoints[ip1L]);
    Point rightPoint = 0.25 * (samplePoints[ip0] + samplePoints[ip1] + samplePoints[ip0R] + samplePoints[ip1R]);

    Vector cx1 = GetDelta(leftPoint, rightPoint, projection);
    // Rotate delta by pi/2
    cx1 = Vector(-cx1.y(), cx1.x());

    Vector cxL = GetDelta(samplePoints[ip0], samplePoints[ip1], projection);
    // Rotate delta by pi/2
    cxL = Vector(-cxL.y(), cxL.x());

    Vector cx0 = -cx1;
    Vector cxR = -cxL;

    double darea = 0.5 * (dot(cx0, x0) + dot(cx1, x1) + dot(cxL, leftPoint) + dot(cxR, rightPoint));

    // !     gradx and grady can be composed

    if (zss(ip0) != constants::missing::doubleValue &&
        zss(ip1) != constants::missing::doubleValue &&
        zss(ip0L) != constants::missing::doubleValue &&
        zss(ip0R) != constants::missing::doubleValue &&
        zss(ip1L) != constants::missing::doubleValue &&
        zss(ip1R) != constants::missing::doubleValue)
    {
        double zL = 0.25 * (zss(ip0) + zss(ip1) + zss(ip0L) + zss(ip1L));
        double zR = 0.25 * (zss(ip0) + zss(ip1) + zss(ip0R) + zss(ip1R));
        gradient[0] = (cx1.x() * z1 + cxL.x() * zL + cx0.x() * z0 + cxR.x() * zR) / darea;
        gradient[1] = (cx1.y() * z1 + cxL.y() * zL + cx0.y() * z0 + cxR.y() * zR) / darea;
    }

    S = 2.0 * cx1;
    dareaL = 0.5 * std::abs(OuterProductTwoSegments(x0, rightPoint, x0, leftPoint, projection));
    dareaR = 0.5 * std::abs(OuterProductTwoSegments(x1, rightPoint, x1, leftPoint, projection));
}

void meshkernel::HessianCalculator::ComputeSampleGradient(const std::vector<Sample>& samplePoints,
                                                          const Projection projection,
                                                          const Hessian& hessian,
                                                          const UInt direction,
                                                          const UInt i,
                                                          const UInt j,
                                                          meshkernel::Vector& gradient,
                                                          meshkernel::Vector& sn,
                                                          double& dareaL,
                                                          double& dareaR)
{

    // const HessianDimension& dimension = hessian.size();

    gradient[0] = 0.0;
    gradient[1] = 0.0;
    sn[0] = 0.0;
    sn[1] = 0.0;

    dareaL = 0.0;
    dareaR = 0.0;

    if (direction == 0)
    {
        //      i-edge gradient at (i+1/2,j) location
        //         control volume:
        //
        //                   L:(i+1/2,j+1/2)
        //                  / \                  |
        //                 /   \                 |
        //          0:(i,j)-----1:(i+1,j)        |
        //                 \   /                 |
        //                  \ /                  |
        //                   R:(i+1/2,j-1/2)

        // UInt ip0 = hessian.get1DIndex(i, j);          // ! pointer to (i,j)
        // UInt ip1 = hessian.get1DIndex(i + 1, j);      // ! pointer to (i+1,j)
        // UInt ip0L = hessian.get1DIndex(i, j + 1);     // ! pointer to (i,j+1)
        // UInt ip0R = hessian.get1DIndex(i, j - 1);     // ! pointer to (i,j-1)
        // UInt ip1L = hessian.get1DIndex(i + 1, j + 1); // ! pointer to (i+1,j+1)
        // UInt ip1R = hessian.get1DIndex(i + 1, j - 1); // ! pointer to (i+1,j-1)

        const auto& dimension = hessian.size();

        UInt ip0 = i + dimension[1] * j;                                  // ! pointer to (i,j)
        UInt ip1 = i + 1 + dimension[1] * j;                              // ! pointer to (i+1,j)
        UInt ip0L = i + dimension[1] * std::min(j + 1, dimension[2]);     // ! pointer to (i,j+1)
        UInt ip0R = i + dimension[1] * std::max(j - 1, 0U);               // ! pointer to (i,j-1)
        UInt ip1L = i + 1 + dimension[1] * std::min(j + 1, dimension[2]); // ! pointer to (i+1,j+1)
        UInt ip1R = i + 1 + dimension[1] * std::max(j - 1, 0U);           // ! pointer to (i+1,j-1)
        ComputeGradient(samplePoints, projection, hessian, ip0, ip1, ip0L, ip0R, ip1L, ip1R, gradient, sn, dareaL, dareaR);
    }
    else if (direction == 1)
    {
        //      j-edge gradient at (i,j+1/2) location
        //        control volume:
        //
        //                   1:(i,j+1)
        //                  / \                  |
        //                 /   \                 |
        //  L:(i-1/2,j+1/2)-----R:(i+1/2,j+1/2)  |
        //                 \   /                 |
        //                  \ /                  |
        //                   0:(i,j)

        // UInt ip0 = hessian.get1DIndex(i, j);
        // UInt ip1 = hessian.get1DIndex(i, j + 1);
        // UInt ip0L = hessian.get1DIndex(i - 1, j);     //              ! pointer to (i-1,j)
        // UInt ip0R = hessian.get1DIndex(i + 1, j);     //              ! pointer to (i+1,j)
        // UInt ip1L = hessian.get1DIndex(i - 1, j + 1); //              ! pointer to (i-1,j+1)
        // UInt ip1R = hessian.get1DIndex(i + 1, j + 1); //              ! pointer to (i+1,j+1)

        const auto& dimension = hessian.size();

        UInt ip0 = i + dimension[1] * j;                                    //              ! pointer to (i,j)
        UInt ip1 = i + dimension[1] * (j + 1);                              //              ! pointer to (i,j+1)
        UInt ip0L = std::max(i - 1, 0U) + dimension[1] * j;                 //              ! pointer to (i-1,j)
        UInt ip0R = std::min(i + 1, dimension[2]) + dimension[1] * j;       //              ! pointer to (i+1,j)
        UInt ip1L = std::max(i - 1, 0U) + dimension[1] * (j + 1);           //              ! pointer to (i-1,j+1)
        UInt ip1R = std::min(i + 1, dimension[2]) + dimension[1] * (j + 1); //              ! pointer to (i+1,j+1)
        ComputeGradient(samplePoints, projection, hessian, ip0, ip1, ip0L, ip0R, ip1L, ip1R, gradient, sn, dareaL, dareaR);
    }
}

void meshkernel::HessianCalculator::ComputeHessian(const std::vector<Sample>& samplePoints,
                                                   const Projection projection,
                                                   Hessian& hessian)
{

    if (hessian.size(1) < 3 || hessian.size(2) < 3)
    {
        return;
    }

    meshkernel::Vector gradientiL;
    meshkernel::Vector gradientiR;
    meshkernel::Vector gradientjL;
    meshkernel::Vector gradientjR;

    meshkernel::Vector SniL;
    meshkernel::Vector SniR;
    meshkernel::Vector SnjL;
    meshkernel::Vector SnjR;

    meshkernel::Vector S;

    Eigen::Matrix2d VV;

    for (UInt i = 1; i < hessian.size(1) - 1; ++i)
    {
        // double af = static_cast<double>(i - 1) / static_cast<double>(std::max(hessian.size(1) - 3, 1));

        for (UInt j = 1; j < hessian.size(2) - 1; ++j)
        {
#if 0
            UInt ip = i + (j - 1) * hessian.size(1);

            // TODO is this some debugging left over?
            if ((std::abs(samplePoints[ip].x - 87270.0) < 1.0e-8) && (std::abs(samplePoints[ip].y - 415570.0) < 1.0e-8))
            {
                continue;
            }
#endif

            double dareaiL = 0.0;
            double dareaiR = 0.0;
            double dareajL = 0.0;
            double dareajR = 0.0;
            double dum; // unused value

            ComputeSampleGradient(samplePoints, projection, hessian, 0, i, j, gradientiR, SniR, dareaiR, dum);

            if (IsEqual(gradientiR[0], constants::missing::doubleValue))
            {
                continue;
            }

            ComputeSampleGradient(samplePoints, projection, hessian, 0, i - 1, j, gradientiL, SniL, dum, dareaiL);

            if (IsEqual(gradientiL[0], constants::missing::doubleValue))
            {
                continue;
            }

            ComputeSampleGradient(samplePoints, projection, hessian, 1, i, j, gradientjR, SnjR, dareajR, dum);

            if (IsEqual(gradientjR[0], constants::missing::doubleValue))
            {
                continue;
            }

            ComputeSampleGradient(samplePoints, projection, hessian, 1, i, j - 1, gradientjL, SnjL, dum, dareajL);

            if (IsEqual(gradientjL[0], constants::missing::doubleValue))
            {
                continue;
            }

            const double area = dareaiL + dareaiR + dareajL + dareajR;
            const double areaInv = 1.0 / area;

            const double zx = (0.5 * (hessian(0, i + 1, j) + hessian(0, i, j)) * SniR[0] - 0.5 * (hessian(0, i - 1, j) + hessian(0, i, j)) * SniL[0] +
                               0.5 * (hessian(0, i, j + 1) + hessian(0, i, j)) * SnjR[0] - 0.5 * (hessian(0, i, j - 1) + hessian(0, i, j)) * SnjL[0]) *
                              areaInv;
            const double zy = (0.5 * (hessian(0, i + 1, j) + hessian(0, i, j)) * SniR[1] - 0.5 * (hessian(0, i - 1, j) + hessian(0, i, j)) * SniL[1] +
                               0.5 * (hessian(0, i, j + 1) + hessian(0, i, j)) * SnjR[1] - 0.5 * (hessian(0, i, j - 1) + hessian(0, i, j)) * SnjL[1]) *
                              areaInv;

            VV(0, 0) = (gradientiR[0] * SniR[0] - gradientiL[0] * SniL[0] + gradientjR[0] * SnjR[0] - gradientjL[0] * SnjL[0]) * areaInv;
            VV(0, 1) = (gradientiR[0] * SniR[1] - gradientiL[0] * SniL[1] + gradientjR[0] * SnjR[1] - gradientjL[0] * SnjL[1]) * areaInv;
            VV(1, 0) = (gradientiR[1] * SniR[0] - gradientiL[1] * SniL[0] + gradientjR[1] * SnjR[0] - gradientjL[1] * SnjL[0]) * areaInv;
            VV(1, 1) = (gradientiR[1] * SniR[1] - gradientiL[1] * SniL[1] + gradientjR[1] * SnjR[1] - gradientjL[1] * SnjL[1]) * areaInv;

            // Eigendecompostion
            Eigen::EigenSolver<Eigen::Matrix2d> eigensolver(VV);

            const Eigen::EigenSolver<Eigen::Matrix2d>::EigenvectorsType& eigenvectors = eigensolver.eigenvectors();
            Eigen::EigenSolver<Eigen::Matrix2d>::EigenvalueType eigenvalues = eigensolver.eigenvalues();

            // Eigen::JacobiSVD svd(VV, Eigen::ComputeFullU | Eigen::ComputeFullV);
            // auto V = svd.matrixV();
            // auto Sigma = svd.singularValues();

            // const auto eigenvalues = Sigma.array().square();
            // const auto eigenvectors = svd.matrixU();

            UInt k = std::abs(eigenvalues[0].real()) > std::abs(eigenvalues[1].real()) ? 0U : 1u;

            hessian(1, i, j) = eigenvectors(0, k).real();
            hessian(2, i, j) = eigenvectors(1, k).real();
            hessian(3, i, j) = eigenvalues[k].real() * area;
            const double uu1 = eigenvectors(0, k).real();
            const double uu2 = eigenvectors(1, k).real();
            const auto firstVectorFirst = eigenvectors(0, 0).real();
            const auto firstVectorSecond = eigenvectors(0, 1).real();
            const auto secondVectorFirst = eigenvectors(k, 0).real();
            const auto secondVectorSecond = eigenvectors(1, 1).real();

            hessian(4, i, j) = -(eigenvectors(k, 0).real() * zx - eigenvectors(k, 1).real() * zy) / (eigenvalues[k].real() + 1.0e-8);
            const double hessianValue = hessian(4, i, j);
            std::cout << hessian(0, i + 1, j) << "  " << hessian(0, i, j) << "  " << hessian(0, i - 1, j) << "  " << hessian(0, i, j) << "  " << hessian(4, i, j) << std::endl;
        }
    }
}

void meshkernel::HessianCalculator::PrepareSampleForHessian(const std::vector<Sample>& samplePoints,
                                                            const Projection projection,
                                                            Hessian& hessian)
{

    // TODO pass in Compute function
    UInt numberOfSmoothingIterations = 0;

    SmoothSamples(samplePoints, numberOfSmoothingIterations, hessian);
    ComputeHessian(samplePoints, projection, hessian);
}

void meshkernel::HessianCalculator::Compute(const std::vector<Sample>& rawSamplePoints,
                                            const Projection projection,
                                            const UInt numX,
                                            const UInt numY,
                                            Hessian& hessian)
{
    std::vector<Sample> samplePoints(rawSamplePoints);

    hessian.resize(5, numX, numY);
    PrepareSampleForHessian(samplePoints, projection, hessian);
}

void meshkernel::HessianCalculator::Compute(const std::vector<Sample>& rawSamplePoints,
                                            const Projection projection,
                                            const UInt numX,
                                            const UInt numY,
                                            std::vector<Sample>& hessianSamples)
{
    hessianSamples = rawSamplePoints;

    Hessian hessian(5, numY, numX);
    PrepareSampleForHessian(hessianSamples, projection, hessian);

    UInt count = 0;

    for (UInt i = 0; i < numX; ++i)
    {
        for (UInt j = 0; j < numY; ++j)
        {

            hessianSamples[count].value = hessian(4, j, i);
            ++count;
        }
    }
}
