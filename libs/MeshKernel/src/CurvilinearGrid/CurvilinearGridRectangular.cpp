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

// include boost
#define BOOST_ALLOW_DEPRECATED_HEADERS
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/formulas/thomas_direct.hpp>
#undef BOOST_ALLOW_DEPRECATED_HEADERS


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
            return std::make_unique<CurvilinearGrid>(ComputeSphericalMercator(numColumns,
            // return std::make_unique<CurvilinearGrid>(ComputeSphericalBoostGrid(numColumns,
            // return std::make_unique<CurvilinearGrid>(ComputeSphericalRgfGrid(numColumns,
            // return std::make_unique<CurvilinearGrid>(ComputeSphericalOnExtension(numColumns,
            // return std::make_unique<CurvilinearGrid>(ComputeSpherical(numColumns,
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
                                                                        const double blockSizeY) const
    {

        lin_alg::Matrix<Point> result = ComputeSphericalOnExtension (numColumns,
                                                                     numRows,
                                                                     originX,
                                                                     originY,
                                                                     angle,
                                                                     blockSizeX,
                                                                     blockSizeY);

        // lin_alg::Matrix<Point> result = ComputeCartesian(numColumns,
        //                                                  numRows,
        //                                                  originX,
        //                                                  originY,
        //                                                  angle,
        //                                                  blockSizeX,
        //                                                  blockSizeY);

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

    Point CurvilinearGridRectangular::RotateByAngle(const double originX, const double originY,
                                                    const double upperRightX, const double upperRightY,
                                                    const double cosAngle,
                                                    const double sinAngle) const
    {

        if (m_projection == Projection::cartesian)
        {
            Point translated(upperRightX - originX, upperRightY - originY);
            return {cosAngle * translated.x - sinAngle * translated.y, sinAngle * translated.x + cosAngle * translated.y};
        }
        else
        {
            Cartesian3DPoint rotationPoint = SphericalToCartesian3D({originX, originY});
            Cartesian3DPoint point3d = SphericalToCartesian3D({upperRightX, upperRightY});

            // Normalize the rotation axis (Rodrigues' formula requires a unit vector)
            // Points are on Earth's surface, so the magnitude is earth_radius
            Cartesian3DPoint k = {rotationPoint.x / constants::geometric::earth_radius,
                                  rotationPoint.y / constants::geometric::earth_radius,
                                  rotationPoint.z / constants::geometric::earth_radius};

            // k \cdot v
            double dot = k.x * point3d.x + k.y * point3d.y + k.z * point3d.z;

            // k \cross v
            Cartesian3DPoint crossProd = VectorProduct(k, point3d);

            // Rodrigues formula: v_rot = v·cos(θ) + (k × v)·sin(θ) + k·(k·v)·(1 - cos(θ))
            Cartesian3DPoint rotatedPoint3d = {point3d.x * cosAngle + crossProd.x * sinAngle + k.x * dot * (1.0 - cosAngle),
                                               point3d.y * cosAngle + crossProd.y * sinAngle + k.y * dot * (1.0 - cosAngle),
                                               point3d.z * cosAngle + crossProd.z * sinAngle + k.z * dot * (1.0 - cosAngle)};

            Point pointOnSphere = Cartesian3DToSpherical(rotatedPoint3d, upperRightX);

            return pointOnSphere;
        }
    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSphericalOnExtension(const int numColumns,
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

        const double cosAngle = std::cos(angle * constants::conversion::degToRad);
        const double sinAngle = std::sin(angle * constants::conversion::degToRad);

        for (Eigen::Index n = 0; n < numN; ++n)
        {

            for (Eigen::Index m = 0; m < numM; ++m)
            {
                result(n, m) = RotateByAngle(originX, originY, result(n, m).x, result(n, m).y, cosAngle, sinAngle);
            }
        }

        return result;
    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSphericalRgfGrid(const int numColumns,
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
                                                         angle,
                                                         blockSizeX,
                                                         blockSizeY);

        // const auto numM = result.cols();
        // const auto numN = result.rows();

        // const double cosAngle = std::cos(angle * constants::conversion::degToRad);
        // const double sinAngle = std::sin(angle * constants::conversion::degToRad);

        // for (Eigen::Index n = 0; n < numN; ++n)
        // {

        //     for (Eigen::Index m = 0; m < numM; ++m)
        //     {
        //         result(n, m) = RotateByAngle(originX, originY, result(n, m).x, result(n, m).y, cosAngle, sinAngle);
        //     }
        // }

        return result;
    }


    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSphericalBoostGrid(const int numColumns,
                                                                                 const int numRows,
                                                                                 const double originX,
                                                                                 const double originY,
                                                                                 const double angle,
                                                                                 const double blockSizeXDeg,
                                                                                 const double blockSizeYDeg) const
    {
        namespace bg = boost::geometry;

        const int numM = numColumns + 1;
        const int numN = numRows + 1;

        using Point2D = bg::model::point<double, 2, bg::cs::geographic<bg::degree>>;

        lin_alg::Matrix<Point> result(numN, numM);

        double angle_y_rad = angle * constants::conversion::degToRad;
        double angle_x_rad = (angle + 90.0) * constants::conversion::degToRad;

        bg::srs::spheroid<double> earth_spheroid(constants::geometric::earth_radius, constants::geometric::earth_radius);
        using thomas_type = bg::formula::thomas_direct<double, true, true, false, false, false>;

        double blockSizeX = blockSizeXDeg * 111000.0 * std::cos (originY * constants::conversion::degToRad);
        double blockSizeY = blockSizeYDeg * 111000.0 * std::cos (originY * constants::conversion::degToRad);

        // using FormulaStrategy = bg::strategy::formula::thomas<double, true, tr, false, false, false>;

        for (Eigen::Index n = 0; n < numN; ++n)
        {

            Point2D current_row_origin;

            if (n == 0) {
                current_row_origin = Point2D(originX, originY);
            } else {
                // Step the next row along the rotated Y-axis vector (heading = angle_y_rad)
                double prev_row_lon_rad = result(n-1, 0).x * bg::math::d2r<double>();
                double prev_row_lat_rad = result(n-1, 0).y * bg::math::d2r<double>();

                auto dir_y = thomas_type::apply (prev_row_lon_rad, prev_row_lat_rad, blockSizeY, angle_y_rad, earth_spheroid);

                current_row_origin = Point2D(dir_y.lon2 * bg::math::r2d<double>(),
                                             dir_y.lat2 * bg::math::r2d<double>());
            }

            for (Eigen::Index m = 0; m < numM; ++m)
            {
                if (m == 0) {
                    result(n, m).x = bg::get<0>(current_row_origin);
                    result(n, m).y = bg::get<1>(current_row_origin);
                } else {

                    [[maybe_unused]] Point cur = result(n, m-1);

                    // Step columns outward along the perpendicular rotated X-axis vector (heading = angle_x_rad)
                    double current_lon_rad = result(n, m-1).x * bg::math::d2r<double>();
                    double current_lat_rad = result(n, m-1).y * bg::math::d2r<double>();


                    auto dir_x = thomas_type::apply (current_lon_rad, current_lat_rad, blockSizeX, angle_x_rad, earth_spheroid);

                    result(n, m).x = dir_x.lon2 * bg::math::r2d<double>();
                    result(n, m).y = dir_x.lat2 * bg::math::r2d<double>();
                }
            }
        }

        return result;
    }

    Cartesian3DPoint CurvilinearGridRectangular::RotateVectorRodrigues(const Cartesian3DPoint& v, const Cartesian3DPoint& k, double theta_rad)  {

        double cos_t = std::cos(theta_rad);
        double sin_t = std::sin(theta_rad);

        Cartesian3DPoint cross = { k.y * v.z - k.z * v.y, k.z * v.x - k.x * v.z, k.x * v.y - k.y * v.x };

        double dot = k.x * v.x + k.y * v.y + k.z * v.z;

        return {
            v.x * cos_t + cross.x * sin_t + k.x * dot * (1.0 - cos_t),
            v.y * cos_t + cross.y * sin_t + k.y * dot * (1.0 - cos_t),
            v.z * cos_t + cross.z * sin_t + k.z * dot * (1.0 - cos_t)
        };
    }

    lin_alg::Matrix<Point> CurvilinearGridRectangular::ComputeSphericalMercator(const int numColumns,
                                                                                const int numRows,
                                                                                const double origin_lon,
                                                                                const double origin_lat,
                                                                                const double rotation_deg,
                                                                                const double d_lon,
                                                                                const double d_lat) const
    {
        const int ny = numColumns + 1;
        const int nx = numRows + 1;

        lin_alg::Matrix<Point> result(ny, nx);

        const double theta_rad = rotation_deg * (M_PI / 180.0);
        Cartesian3DPoint k = SphericalToCartesian3D(Point(origin_lon, origin_lat));

        double kLength = std::sqrt (k.x * k.x + k.y * k.y + k.z * k.z);

        k.x /= kLength;
        k.y /= kLength;
        k.z /= kLength;

        // 1. Transform starting origin latitude into Conformal Mercator space
        double origin_lat_rad = origin_lat * (M_PI / 180.0);
        double origin_y_mercator = std::log(std::tan(M_PI / 4.0 + origin_lat_rad / 2.0));

        // 2. Compute the precise scaling factor for your rectangular dimensions.
        // Instead of forcing square steps, we scale the vertical Mercator increment
        // by the exact ratio of your desired d_lat to d_lon.
        double d_step_lon_rad = d_lon * (M_PI / 180.0);
        double d_step_lat_conformal = d_step_lon_rad * (d_lat / d_lon);

        for (int j = 0; j < ny; ++j) {
            // Step along the scaled conformal Y axis
            double current_y_mercator = origin_y_mercator + (j * d_step_lat_conformal);

            // Inverse Mercator equation back to physical latitude
            double unrotated_lat_rad = 2.0 * std::atan(std::exp(current_y_mercator)) - M_PI / 2.0;
            double unrotated_lat = unrotated_lat_rad * (180.0 / M_PI);

            for (int i = 0; i < nx; ++i) {
                // Step uniformly along the X axis
                double unrotated_lon = origin_lon + (i * d_lon);

                // Convert conformal unrotated point into 3D space
                Cartesian3DPoint p_initial = SphericalToCartesian3D(Point(unrotated_lon, unrotated_lat));

                // Apply 3D Rodrigues rotation
                Cartesian3DPoint p_rotated = RotateVectorRodrigues(p_initial, k, theta_rad);

                // Project back to degrees and update the Eigen Matrix
                double final_lon, final_lat;

                final_lat = std::asin(std::max(-1.0, std::min (1.0, p_rotated.z))) * (180.0 / M_PI);
                final_lon = std::atan2(p_rotated.y, p_rotated.x) * (180.0 / M_PI);

                // std::cout << current_y_mercator << "  "<< k.x << ", " << k.y << ", " << k.z << "   " << final_lat << "  " << final_lon << "  " << rotation_deg << "  " << d_lon << "   " << d_lat << std::endl;
                // std::cout << current_y_mercator << "  "<< p_initial.x << ", " << p_initial.y << ", " << p_initial.z << "   " << final_lat << "  " << final_lon << "  " << rotation_deg << "  " << d_lon << "   " << d_lat << std::endl;
                // std::cout << current_y_mercator << "  "<< p_rotated.x << ", " << p_rotated.y << ", " << p_rotated.z << "   " << final_lat << "  " << final_lon << "  " << rotation_deg << "  " << d_lon << "   " << d_lat << "   " << theta_rad << std::endl;

                //CartesianToLatLon(p_rotated, final_lon, final_lat);

                result(j, i).x = final_lon;
                result(j, i).y = final_lat;
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

        const double cosAngle = std::cos(-angle * constants::conversion::degToRad);
        const double sinAngle = std::sin(-angle * constants::conversion::degToRad);

        // rotate the upper right, by -angle so that the grid is aligned with the axis
        Point rotatedUpperRight = RotateByAngle(originX, originY, upperRightX, upperRightY, cosAngle, sinAngle);

        // Now the number of cells in each direction can be computed.
        const int numColumns = static_cast<int>(std::ceil((rotatedUpperRight.x - originX) / blockSizeX));
        const int numRows = ComputeNumRows(originY, rotatedUpperRight.y, blockSizeX, blockSizeY, m_projection);

        if (numColumns <= 0)
        {
            throw AlgorithmError("Number of columns cannot be <= 0");
        }

        if (numRows <= 0)
        {
            throw AlgorithmError("Number of rows cannot be <= 0");
        }

        if (m_projection == Projection::spherical)
        {
            auto grid = std::make_unique<CurvilinearGrid>(ComputeSphericalOnExtension(numColumns,
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
