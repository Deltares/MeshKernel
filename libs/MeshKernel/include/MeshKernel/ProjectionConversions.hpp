//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2023.
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

#include <string>

#include <boost/geometry/core/coordinate_system.hpp>
#include <boost/geometry/srs/epsg.hpp>

#include "MeshKernel/Mesh.hpp"
#include "MeshKernel/Vector.hpp"

namespace meshkernel
{

    /// @brief Namespace alias for boost::geometry
    namespace bg = boost::geometry;

    /// @brief Converts points from spherical to Cartesian coordinate system.
    template <typename ProjectionConversion>
    class ConvertSphericalToCartesianBase
    {
    public:
        /// @brief point in longitude-latitude space
        using LongLat = bg::model::d2::point_xy<double, bg::cs::geographic<bg::degree>>;

        /// @brief Point in x-y space
        using UTM = bg::model::d2::point_xy<double, Projection>;

        /// @brief Constructor with projection
        ConvertSphericalToCartesianBase(const ProjectionConversion& proj) : m_projection(proj) {}

        /// @brief Default destructor
        virtual ~ConvertSphericalToCartesianBase() = default;

        /// @brief The coordinate system of the point parameter to the conversion operation.
        Projection SourceProjection() const
        {
            return Projection::spherical;
        }

        /// @brief The coordinate system of the point result of the conversion operation.
        Projection TargetProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief Apply the conversion of a point in Spherical coordinate system to Cartesian
        Point operator()(const Point& pnt) const
        {
            LongLat longLat(pnt.x, pnt.y);
            UTM utm{0.0, 0.0};
            m_projection.forward(longLat, utm);

            Point result(utm.x(), utm.y());
            return result;
        }

    private:
        /// @brief The projection conversion object.
        ProjectionConversion m_projection;
    };

    /// @brief Converts points from spherical to Cartesian coordinate system using an ESPG code.
    template <const int EpsgCode>
    class ConvertSphericalToCartesianEPSG : public ConvertSphericalToCartesianBase<bg::srs::projection<bg::srs::static_epsg<EpsgCode>>>
    {
    public:
        /// @brief The EPSG projection
        using EpsgProjection = bg::srs::projection<bg::srs::static_epsg<EpsgCode>>;

        /// @brief Construct spherical to Cartesian with an EPSG code
        ConvertSphericalToCartesianEPSG() : ConvertSphericalToCartesianBase<bg::srs::projection<bg::srs::static_epsg<EpsgCode>>>(EpsgProjection()) {}
    };

    /// @brief Converts points from spherical to Cartesian coordinate system.
    class ConvertSphericalToCartesian : public ConvertSphericalToCartesianBase<bg::srs::projection<>>
    {
    public:
        /// @brief Construct spherical to Cartesian with an zone string
        ConvertSphericalToCartesian(const std::string& zone) : ConvertSphericalToCartesianBase<bg::srs::projection<>>(bg::srs::proj4(zone)) {}
    };

    //--------------------------------

    /// @brief Converts points from spherical to Cartesian coordinate system.
    template <typename ProjectionConversion>
    class ConvertCartesianToSphericalBase
    {
    public:
        /// @brief point in longitude-latitude space
        using LongLat = bg::model::d2::point_xy<double, bg::cs::geographic<bg::degree>>;

        /// @brief Point in x-y space
        using UTM = bg::model::d2::point_xy<double>;

        /// @brief Constructor with projection
        ConvertCartesianToSphericalBase(const ProjectionConversion& proj) : m_projection(proj) {}

        /// @brief Default destructor
        virtual ~ConvertCartesianToSphericalBase() = default;

        /// @brief The coordinate system of the point parameter to the conversion operation.
        Projection SourceProjection() const
        {
            return Projection::cartesian;
        }

        /// @brief The coordinate system of the point result of the conversion operation.
        Projection TargetProjection() const
        {
            return Projection::spherical;
        }

        /// @brief Apply the conversion of a point in Cartesian coordinate system to spherical
        Point operator()(const Point& pnt) const
        {
            UTM utm{pnt.x, pnt.y};
            LongLat longLat{0.0, 0.0};
            m_projection.inverse(utm, longLat);

            Point result(longLat.x(), longLat.y());
            return result;
        }

    private:
        /// @brief The projection conversion object.
        ProjectionConversion m_projection;
    };

    /// @brief Converts points from spherical to Cartesian coordinate system using an EPSG code.
    template <const int EpsgCode>
    class ConvertCartesianToSphericalEPSG : public ConvertCartesianToSphericalBase<boost::geometry::srs::projection<boost::geometry::srs::static_epsg<EpsgCode>>>
    {
    public:
        /// @brief The EPSG projection
        using EpsgProjection = boost::geometry::srs::projection<boost::geometry::srs::static_epsg<EpsgCode>>;

        /// @brief Construct spherical to Cartesian with an EPSG code
        ConvertCartesianToSphericalEPSG() : ConvertCartesianToSphericalBase<boost::geometry::srs::projection<boost::geometry::srs::static_epsg<EpsgCode>>>(EpsgProjection()) {}
    };

    /// @brief Converts points from spherical to Cartesian coordinate system.
    class ConvertCartesianToSpherical : public ConvertCartesianToSphericalBase<bg::srs::projection<>>
    {
    public:
        /// @brief Construct spherical to Cartesian with an zone string
        ConvertCartesianToSpherical(const std::string& zone) : ConvertCartesianToSphericalBase<bg::srs::projection<>>(bg::srs::proj4(zone)) {}
    };

} // namespace meshkernel
