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

#pragma once

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>

#include <stdexcept>
#include <vector>

namespace meshkernelapi
{

    /// @brief Converts a GeometryList to a vector<Point>
    /// @param[in] geometryListIn The geometry input list to convert
    /// @returns The converted vector of points
    static std::vector<meshkernel::Point> ConvertGeometryListToPointVector(const GeometryList& geometryListIn)
    {
        std::vector<meshkernel::Point> result;
        if (geometryListIn.num_coordinates == 0)
        {
            return result;
        }
        result.reserve(geometryListIn.num_coordinates);
        result.emplace_back(geometryListIn.coordinates_x[0], geometryListIn.coordinates_y[0]);
        // remove consecutive duplicated point leading to 0 edge length
        for (auto i = 1; i < geometryListIn.num_coordinates; ++i)
        {
            if (meshkernel::IsEqual(geometryListIn.coordinates_x[i], result.back().x) &&
                meshkernel::IsEqual(geometryListIn.coordinates_y[i], result.back().y))
            {
                continue;
            }
            result.emplace_back(geometryListIn.coordinates_x[i], geometryListIn.coordinates_y[i]);
        }
        return result;
    }

    /// @brief Converts a GeometryList to a vector<vector<Point>>
    /// @param[in] geometryListIn The geometry input list to convert
    /// @returns The converted vector of points
    /// @return The resulting vector of vector points
    static std::vector<std::vector<meshkernel::Point>> ConvertGeometryListToVectorOfPointVectors(GeometryList const& geometryListIn)
    {
        std::vector<std::vector<meshkernel::Point>> result;
        std::vector<meshkernel::Point> chunk;
        chunk.reserve(geometryListIn.num_coordinates);
        for (auto i = 0; i < geometryListIn.num_coordinates; ++i)
        {

            if (meshkernel::IsEqual(geometryListIn.coordinates_x[i], geometryListIn.geometry_separator))
            {
                result.emplace_back(chunk);
                chunk.clear();
            }
            else
            {
                chunk.emplace_back(geometryListIn.coordinates_x[i], geometryListIn.coordinates_y[i]);
            }
        }

        if (!chunk.empty())
        {
            result.emplace_back(chunk);
        }

        return result;
    }

    /// @brief Converts an a vector<T> separated by \p separator to a vector<vector<T>>
    /// @param[in]  input The data of the input array
    /// @param[in]  separator  The separator value
    /// @return The resulting vector of vectors
    template <typename T>
    static std::vector<std::vector<T>> ConvertVectorToVectorOfVectors(const std::vector<T> input, T separator)
    {
        std::vector<std::vector<T>> result;
        std::vector<T> chunk;
        for (size_t i = 0; i < input.size(); ++i)
        {
            if (meshkernel::IsEqual(input[i], separator))
            {
                result.emplace_back(chunk);
                chunk.clear();
            }
            else
            {
                chunk.emplace_back(input[i]);
            }
        }
        if (!chunk.empty())
        {
            result.emplace_back(chunk);
        }
        return result;
    }

    /// @brief Converts an int array to a vector<bool>
    /// @param[in]  inputArray The data of the input array
    /// @param[in]  inputSize  The size of the input array
    /// @return The resulting vector of bool
    static std::vector<bool> ConvertIntegerArrayToBoolVector(const int inputArray[], size_t inputSize)
    {
        std::vector<bool> result(inputSize);
        for (size_t i = 0; i < inputSize; ++i)
        {
            switch (inputArray[i])
            {
            case 0:
                result[i] = false;
                break;
            case 1:
                result[i] = true;
                break;
            default:
                throw std::invalid_argument("MeshKernel: Invalid 1D mask.");
            }
        }
        return result;
    }

    /// @brief Converts a GeometryList to a vector<Sample>
    /// @param[in] geometryListIn The geometry input list to convert
    /// @return The converted vector of samples
    static std::vector<meshkernel::Sample> ConvertGeometryListToSampleVector(const GeometryList& geometryListIn)
    {
        if (geometryListIn.num_coordinates == 0)
        {
            throw std::invalid_argument("MeshKernel: The samples are empty.");
        }
        std::vector<meshkernel::Sample> result;
        result.reserve(geometryListIn.num_coordinates);

        for (auto i = 0; i < geometryListIn.num_coordinates; ++i)
        {
            result.push_back({geometryListIn.coordinates_x[i], geometryListIn.coordinates_y[i], geometryListIn.values[i]});
        }
        return result;
    }

    /// @brief Converts a vector<Point> to a GeometryList
    /// @param[in]  pointVector The point vector to convert
    /// @param[out] result      The converted geometry list
    static void ConvertPointVectorToGeometryList(std::vector<meshkernel::Point> const& pointVector, GeometryList& result)
    {
        if (pointVector.size() < static_cast<size_t>(result.num_coordinates))
        {
            throw std::invalid_argument("MeshKernel: Invalid memory allocation, the point-vector size is smaller than the number of coordinates in the result vector.");
        }

        for (auto i = 0; i < result.num_coordinates; ++i)
        {
            result.coordinates_x[i] = pointVector[i].x;
            result.coordinates_y[i] = pointVector[i].y;
        }
    }

    /// @brief Converts valuesCoordinates and values to a GeometryList result
    /// @param[in]  valuesCoordinates The vector of coordinates
    /// @param[in]  values            The vector of values at the coordinate locations
    /// @param[out] result            The converted geometry list
    static void ConvertSampleVectorToGeometryList(std::vector<meshkernel::Point> const& valuesCoordinates, std::vector<double> const& values, GeometryList& result)
    {
        if (valuesCoordinates.size() != values.size())
        {
            throw std::invalid_argument("MeshKernel: The size of the valuesCoordinates-vector is not equal to the size of the values-vector");
        }

        if (values.size() < static_cast<size_t>(result.num_coordinates))
        {
            throw std::invalid_argument("MeshKernel: Invalid memory allocation, the value-vector size is smaller than the number of coordinates in the result vector.");
        }

        for (auto i = 0; i < result.num_coordinates; ++i)
        {
            result.coordinates_x[i] = valuesCoordinates[i].x;
            result.coordinates_y[i] = valuesCoordinates[i].y;
            result.values[i] = values[i];
        }
    }

    /// @brief Sets splines from a geometry list
    /// @param[in]  geometryListIn The input geometry list
    /// @param[out] spline         The spline which will be set
    static void SetSplines(const GeometryList& geometryListIn, meshkernel::Splines& spline)
    {
        if (geometryListIn.num_coordinates == 0)
        {
            return;
        }

        const auto splineCornerPoints = ConvertGeometryListToPointVector(geometryListIn);

        const auto indices = FindIndices(splineCornerPoints, 0, splineCornerPoints.size(), meshkernel::constants::missing::doubleValue);

        for (const auto& index : indices)
        {
            const auto& [startIndex, endIndex] = index;
            const auto size = endIndex - startIndex + 1;
            if (size > 0)
            {
                spline.AddSpline(splineCornerPoints, startIndex, size);
            }
        }
    }

    /// @brief Sets dimensions members of meshkernelapi::Mesh2D instance
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetMesh2dApiDimensions(std::shared_ptr<meshkernel::Mesh> mesh2d, Mesh2D& mesh2dApi)
    {
        size_t num_face_nodes = 0;
        for (size_t f = 0; f < mesh2d->GetNumFaces(); f++)
        {
            num_face_nodes += mesh2d->m_facesNodes[f].size();
        }

        mesh2dApi.num_face_nodes = static_cast<int>(num_face_nodes);
        mesh2dApi.num_faces = static_cast<int>(mesh2d->GetNumFaces());
        mesh2dApi.num_nodes = static_cast<int>(mesh2d->GetNumNodes());
        mesh2dApi.num_edges = static_cast<int>(mesh2d->GetNumEdges());
    }

    /// @brief Sets the meshkernelapi::Mesh2D data
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetMesh2dApiData(std::shared_ptr<meshkernel::Mesh2D> mesh2d, Mesh2D& mesh2dApi)
    {
        mesh2d->ComputeEdgesCenters();
        for (size_t n = 0; n < mesh2d->GetNumNodes(); n++)
        {
            mesh2dApi.node_x[n] = mesh2d->m_nodes[n].x;
            mesh2dApi.node_y[n] = mesh2d->m_nodes[n].y;
        }

        for (size_t edgeIndex = 0; edgeIndex < mesh2d->GetNumEdges(); edgeIndex++)
        {
            mesh2dApi.edge_x[edgeIndex] = mesh2d->m_edgesCenters[edgeIndex].x;
            mesh2dApi.edge_y[edgeIndex] = mesh2d->m_edgesCenters[edgeIndex].y;
            mesh2dApi.edge_nodes[edgeIndex * 2] = static_cast<int>(mesh2d->m_edges[edgeIndex].first);
            mesh2dApi.edge_nodes[edgeIndex * 2 + 1] = static_cast<int>(mesh2d->m_edges[edgeIndex].second);
        }

        int faceIndex = 0;
        for (size_t f = 0; f < mesh2d->GetNumFaces(); f++)
        {
            mesh2dApi.face_x[f] = mesh2d->m_facesMassCenters[f].x;
            mesh2dApi.face_y[f] = mesh2d->m_facesMassCenters[f].y;
            mesh2dApi.nodes_per_face[f] = static_cast<int>(mesh2d->m_facesNodes[f].size());
            for (size_t n = 0; n < mesh2d->m_facesNodes[f].size(); ++n)
            {
                mesh2dApi.face_nodes[faceIndex] = static_cast<int>(mesh2d->m_facesNodes[f][n]);
                faceIndex++;
            }
        }
    }

    /// @brief Sets the meshkernelapi::CurvilinearGrid data
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetCurvilinearGridApiData(std::shared_ptr<meshkernel::CurvilinearGrid> curvilinearGrid,
                                          CurvilinearGrid& curvilinearGridApi)
    {
        curvilinearGrid->ComputeEdgesCenters();
        for (size_t n = 0; n < curvilinearGrid->GetNumNodes(); n++)
        {
            curvilinearGridApi.node_x[n] = curvilinearGrid->m_nodes[n].x;
            curvilinearGridApi.node_y[n] = curvilinearGrid->m_nodes[n].y;
        }
    }

    /// @brief Sets a meshkernelapi::Mesh1D data
    /// @param[in]  mesh1d           The input meshkernel::Mesh1D instance
    /// @param[out] mesh1dApi        The output meshkernelapi::Mesh1D instance
    static void SetMesh1dApiData(std::shared_ptr<meshkernel::Mesh1D> mesh1d,
                                 Mesh1D& mesh1dApi)
    {
        for (size_t n = 0; n < mesh1d->GetNumNodes(); n++)
        {
            mesh1dApi.node_x[n] = mesh1d->m_nodes[n].x;
            mesh1dApi.node_y[n] = mesh1d->m_nodes[n].y;
        }

        size_t edgeIndex = 0;
        for (size_t e = 0; e < mesh1d->GetNumEdges(); e++)
        {
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d->m_edges[e].first);
            edgeIndex++;
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d->m_edges[e].second);
            edgeIndex++;
        }
    }

} // namespace meshkernelapi
