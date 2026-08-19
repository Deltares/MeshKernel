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

#include "MeshKernel/Exceptions.hpp"

#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRectangular.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Polygons.hpp>
#include <MeshKernel/RangeCheck.hpp>

#include <cmath>

namespace meshkernel
{

    CurvilinearGridRectangular::CurvilinearGridRectangular(Projection projection) : m_projection(projection)
    {
        if (m_projection != Projection::cartesian && m_projection != Projection::spherical)
        {
            throw meshkernel::NotImplementedError("Projection value: {} not supported", static_cast<int>(m_projection));
        }
    }

    std::unique_ptr<CurvilinearGrid> CurvilinearGridRectangular::Compute(const int numColumns,
                                                                         const int numRows,
                                                                         const double originX,
                                                                         const double originY,
                                                                         const double angle,
                                                                         const double blockSizeX,
                                                                         const double blockSizeY) const
    {
        range_check::CheckGreater(numColumns, 0, "Number of columns");
        range_check::CheckGreater(numRows, 0, "Number of rows");
        range_check::CheckInOpenInterval(angle, {-90.0, 90.0}, "Grid angle");
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        if (m_projection == Projection::spherical)
        {
            return std::make_unique<CurvilinearGrid>(ComputeSpherical(numColumns,
                                                                      numRows,
                                                                      originX,
                                                                      originY,
                                                                      angle,
                                                                      blockSizeX,
                                                                      blockSizeY),
                                                     m_projection);
        }
        if (m_projection == Projection::cartesian)
        {
            return std::make_unique<CurvilinearGrid>(ComputeCartesian(numColumns,
                                                                      numRows,
                                                                      originX,
                                                                      originY,
                                                                      angle,
                                                                      blockSizeX,
                                                                      blockSizeY),
                                                     m_projection);
        }
        throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeCartesian(const int numColumns,
                                                                        const int numRows,
                                                                        const double originX,
                                                                        const double originY,
                                                                        const double angle,
                                                                        const double blockSizeX,
                                                                        const double blockSizeY)

    {
        const auto angleInRad = angle * constants::conversion::degToRad;
        const auto cosineAngle = std::cos(angleInRad);
        const auto sinAngle = std::sin(angleInRad);

        const auto numM = numColumns + 1;
        const auto numN = numRows + 1;

        lin_alg::Matrix<Point> result(numN, numM);
        const auto blockSizeXByCos = blockSizeX * cosineAngle;
        const auto blockSizeYbySin = blockSizeY * sinAngle;
        const auto blockSizeXBySin = blockSizeX * sinAngle;
        const auto blockSizeYByCos = blockSizeY * cosineAngle;
        for (Eigen::Index n = 0; n < result.rows(); ++n)
        {
            for (Eigen::Index m = 0; m < result.cols(); ++m)
            {
                const double newPointXCoordinate = originX + m * blockSizeXByCos - n * blockSizeYbySin;
                const double newPointYCoordinate = originY + m * blockSizeXBySin + n * blockSizeYByCos;
                result(n, m) = {newPointXCoordinate, newPointYCoordinate};
            }
        }
        return result;
    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSpherical(const int numColumns,
                                                                        const int numRows,
                                                                        const double originX,
                                                                        const double originY,
                                                                        const double angle,
                                                                        const double blockSizeX,
                                                                        const double blockSizeY)
    {
        lin_alg::Matrix<Point> result = ComputeCartesian(numColumns,
                                                         numRows,
                                                         originX,
                                                         originY,
                                                         angle,
                                                         blockSizeX,
                                                         blockSizeY);

        const auto numM = result.cols();
        const auto numN = result.rows();
        const double latitudePoles = 90.0;
        const double aspectRatio = blockSizeY / blockSizeX;

        bool onPoles = false;
        Eigen::Index lastRowOnPole = numM;

        for (Eigen::Index n = 1; n < numN; ++n)
        {

            for (Eigen::Index m = 0; m < numM; ++m)
            {
                double latitude = ComputeLatitudeIncrementWithAdjustment(blockSizeX, aspectRatio, result(n - 1, m).y);

                result(n, m).y = latitude;

                if (const double latitudeAbs = std::abs(latitude); IsEqual(latitudeAbs, latitudePoles) || latitudeAbs >= latitudePoles)
                {
                    onPoles = true;
                    lastRowOnPole = n;
                }
            }

            if (onPoles)
            {
                if (lastRowOnPole + 1 < result.rows())
                {
                    lin_alg::EraseRows(result, lastRowOnPole + 1, result.rows() - 1);
                }
                break;
            }
        }

        return result;
    }

    double CurvilinearGridRectangular::ComputeLatitudeIncrementWithAdjustment(double blockSize, double aspectRatio, double latitude)
    {

        // When the real distance along the latitude becomes smaller than minimumDistance
        // and the location is close to the poles, snap the next point to the poles.
        const double minimumDistance = 1000.0;

        // The latitude defining close to poles
        const double latitudeCloseToPole = 89.0;

        // The haversine function is defined as:
        //
        // dlon = abs(lon2 - lon1)
        // dlat = abs(lat2 - lat1)
        // a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
        // c = 2 * asin(sqrt(a))
        // dist = radius * c
        //
        // Now, assuming the angle with which the mesh is to be rotated is zero
        // the longitudes lon1 and lon2 are the separated by the blockSizeX
        // the latitudes of the two points are the same, so lat1 - lat2 = 0 => sin (dlat / 2) = 0
        // We use the current latitude to compute the distance to the next point

        double sinLon = std::sin(0.5 * blockSize * constants::conversion::degToRad);
        double cosLat = std::cos(latitude * constants::conversion::degToRad);
        double a = sinLon * sinLon * cosLat * cosLat;
        double c = 2.0 * std::asin(std::sqrt(a));
        c *= aspectRatio;

        double distance = c * constants::geometric::earth_radius;

        double computedLatitude = c * constants::conversion::radToDeg + latitude;

        if (std::abs(computedLatitude) > latitudeCloseToPole && distance < minimumDistance)
        {
            computedLatitude = std::copysign(1.0, computedLatitude) * 90.0;
        }

        return computedLatitude;
    }

    int CurvilinearGridRectangular::ComputeNumRows(double minY,
                                                   double maxY,
                                                   double blockSizeX,
                                                   double blockSizeY,
                                                   Projection projection)
    {
        if (blockSizeY > std::abs(maxY - minY))
        {
            throw AlgorithmError("blockSizeY cannot be larger than mesh height");
        }

        if (projection == Projection::cartesian)
        {
            const int numM = static_cast<int>(std::ceil(std::abs(maxY - minY) / blockSizeY));
            return std::max(numM, 1);
        }

        double currentLatitude = minY;
        int result = 0;
        const double latitudePoles = 90.0;

        double aspectRatio = blockSizeY / blockSizeX;

        while (currentLatitude < maxY)
        {
            currentLatitude = ComputeLatitudeIncrementWithAdjustment(blockSizeX, aspectRatio, currentLatitude);

            result += 1;

            if (IsEqual(std::abs(currentLatitude), latitudePoles))
            {
                break;
            }
        }

        return result;
    }

    std::unique_ptr<CurvilinearGrid> CurvilinearGridRectangular::Compute(const double angle,
                                                                         const double blockSizeX,
                                                                         const double blockSizeY,
                                                                         std::shared_ptr<Polygons> polygons,
                                                                         UInt polygonIndex) const
    {

        range_check::CheckInOpenInterval(angle, {-90.0, 90.0}, "Grid angle");
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        if (polygons->IsEmpty())
        {
            throw AlgorithmError("Enclosures list is empty.");
        }

        if (polygons->GetProjection() != m_projection)
        {
            throw AlgorithmError("Polygon projection ({}) is not equal to curvilinear grid projection ({})",
                                 ProjectionToString(polygons->GetProjection()),
                                 ProjectionToString(m_projection));
        }

        // Compute the bounding box
        const auto boundingBox = polygons->GetBoundingBox(polygonIndex);
        const auto referencePoint = boundingBox.MassCentre();

        // Compute the max size
        const auto maxSize = std::max(boundingBox.Width(), boundingBox.Height());

        // Compute the lower left and upper right corners
        const Point lowerLeft(referencePoint.x - maxSize, referencePoint.y - maxSize);
        const Point upperRight(referencePoint.x + maxSize, referencePoint.y + maxSize);

        // Compute the number of rows and columns
        const int numColumns = std::max(static_cast<int>(std::ceil(std::abs(upperRight.x - lowerLeft.x) / blockSizeX)), 1);
        const int numRows = ComputeNumRows(lowerLeft.y, upperRight.y, blockSizeX, blockSizeY, m_projection);

        // Rotated the lower left corner
        const auto lowerLeftMergedRotated = Rotate(lowerLeft, angle, referencePoint);

        // Set the origin
        const double originX = lowerLeftMergedRotated.x;
        const double originY = lowerLeftMergedRotated.y;

        if (m_projection == Projection::spherical)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeSpherical(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           angle,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);
            grid->Delete(polygons, polygonIndex);
            return grid;
        }
        if (m_projection == Projection::cartesian)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeCartesian(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           angle,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);
            grid->Delete(polygons, polygonIndex);
            return grid;
        }

        throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
    }

    std::unique_ptr<CurvilinearGrid> CurvilinearGridRectangular::Compute(const double originX,
                                                                         const double originY,
                                                                         const double blockSizeX,
                                                                         const double blockSizeY,
                                                                         const double upperRightX,
                                                                         const double upperRightY) const
    {
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        const int numColumns = static_cast<int>(std::ceil((upperRightX - originX) / blockSizeX));
        if (numColumns <= 0)
        {
            throw AlgorithmError("Number of columns cannot be <= 0");
        }

        const int numRows = ComputeNumRows(originY, upperRightY, blockSizeX, blockSizeY, m_projection);

        if (m_projection == Projection::spherical)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeSpherical(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           0.0,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);

            return grid;
        }
        if (m_projection == Projection::cartesian)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeCartesian(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           0.0,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);
            return grid;
        }
        throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
    }

    Point CurvilinearGridRectangular::RotateUpperRightByAngle (const double originX, const double originY, const double upperRightX, const double  upperRightY, const double angle) const
    {

        const double cosA = std::cos(angle * constants::conversion::degToRad);
        const double sinA = std::sin(angle * constants::conversion::degToRad);

        if (m_projection == Projection::cartesian)
        {
            Point translated (upperRightX - originX, upperRightY - originY);
            return {cosA * translated.x - sinA * translated.y, sinA * translated.x + cosA * translated.y};
        }
        else
        {
            Cartesian3DPoint rotationPoint = SphericalToCartesian3D ({originX, originY});
            Cartesian3DPoint point3d = SphericalToCartesian3D ({upperRightX, upperRightY});

            // Normalize the rotation axis (Rodrigues' formula requires a unit vector)
            // Points are on Earth's surface, so the magnitude is earth_radius
            Cartesian3DPoint k = {rotationPoint.x / constants::geometric::earth_radius,
                                  rotationPoint.y / constants::geometric::earth_radius,
                                  rotationPoint.z / constants::geometric::earth_radius};

            // std::cout << "rotationPoint " << rotationPoint.x << ", " << rotationPoint.y << ", " << rotationPoint.z << std::endl;
            // std::cout << "point3d       " << point3d.x << ", " << point3d.y << ", " << point3d.z << std::endl;
            // std::cout << "k (unit axis) " << k.x << ", " << k.y << ", " << k.z << std::endl;

            // k · v
            double dot = k.x * point3d.x + k.y * point3d.y + k.z * point3d.z;

            // k × v
            Cartesian3DPoint crossProd = VectorProduct(k, point3d);
            // std::cout << "crossProd    " << crossProd.x << ", " << crossProd.y << ", " << crossProd.z << std::endl;

            // Rodrigues' formula: v_rot = v·cos(θ) + (k × v)·sin(θ) + k·(k·v)·(1 - cos(θ))
            Cartesian3DPoint rotatedPoint3d = {point3d.x * cosA + crossProd.x * sinA + k.x * dot * (1.0 - cosA),
                                               point3d.y * cosA + crossProd.y * sinA + k.y * dot * (1.0 - cosA),
                                               point3d.z * cosA + crossProd.z * sinA + k.z * dot * (1.0 - cosA)};

            // std::cout << "rotatedPoint3d " << rotatedPoint3d.x << ", " << rotatedPoint3d.y << ", " << rotatedPoint3d.z << std::endl;

            Point pointOnSphere = Cartesian3DToSpherical (rotatedPoint3d, upperRightX);

            std::cout << "pointOnSphere "<< pointOnSphere.x << ",  " << pointOnSphere.y << std::endl;

            return pointOnSphere;
        }

    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSpherical2(const int numColumns,
                                                                        const int numRows,
                                                                        const double originX,
                                                                        const double originY,
                                                                        const double angle,
                                                                        const double blockSizeX,
                                                                        const double blockSizeY) const
    {
        lin_alg::Matrix<Point> result = ComputeCartesian(numColumns,
                                                         numRows,
                                                         originX,
                                                         originY,
                                                         0.0 * angle,
                                                         blockSizeX,
                                                         blockSizeY);

        const auto numM = result.cols();
        const auto numN = result.rows();

        for (Eigen::Index n = 0; n < numN; ++n)
        {

            for (Eigen::Index m = 0; m < numM; ++m)
            {
                result (n,m) = RotateUpperRightByAngle (originX, originY, result (n,m).x, result (n,m).y, angle);
            }
        }


        return result;
    }


    std::unique_ptr<CurvilinearGrid> CurvilinearGridRectangular::Compute(const double originX,
                                                                         const double originY,
                                                                         const double blockSizeX,
                                                                         const double blockSizeY,
                                                                         const double upperRightX,
                                                                         const double upperRightY,
                                                                         const double angle) const
    {
        range_check::CheckGreater(blockSizeX, 0.0, "X block size");
        range_check::CheckGreater(blockSizeY, 0.0, "Y block size");

        Point rotatedUpperRight = RotateUpperRightByAngle (originX, originY, upperRightX, upperRightY, -angle);

        std::cout << "rotated distance: " << rotatedUpperRight.x - originX << ",  " << rotatedUpperRight.y - originY << std::endl;
        std::cout << "rotated deltas  : " << (rotatedUpperRight.x - originX) / blockSizeX << ",  " << (rotatedUpperRight.y - originY) / blockSizeY << std::endl;


        const int numColumns = static_cast<int>(std::ceil((rotatedUpperRight.x - originX) / blockSizeX));

        if (numColumns <= 0)
        {
            throw AlgorithmError("Number of columns cannot be <= 0");
        }

        std::cout << "rotatedUpperRight "<< rotatedUpperRight.x << ", " << rotatedUpperRight.y << std::endl;

        const int numRows = static_cast<int>(std::ceil((rotatedUpperRight.y - originY) / blockSizeY));
        // const int numRows = ComputeNumRows(originY, rotatedUpperRight.y, blockSizeX, blockSizeY, m_projection);
        std::cout << "num rows, cols "<< numRows << "  " << numColumns << "  "<< originX << ", " << originY << " "<< rotatedUpperRight.x << ", " << rotatedUpperRight.y << std::endl;

        if (m_projection == Projection::spherical)
        {
            // double deltax = (rotatedUpperRight.x - originX) / blockSizeX;
            // double bsx =
            auto grid = std::make_unique<CurvilinearGrid>(ComputeSpherical(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           angle,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);

            return grid;
        }
        if (m_projection == Projection::cartesian)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeCartesian(numColumns,
                                                                           numRows,
                                                                           originX,
                                                                           originY,
                                                                           angle,
                                                                           blockSizeX,
                                                                           blockSizeY),
                                                          m_projection);
            return grid;
        }
        throw NotImplementedError("Projection value {} not supported", static_cast<int>(m_projection));
    }

} // namespace meshkernel
