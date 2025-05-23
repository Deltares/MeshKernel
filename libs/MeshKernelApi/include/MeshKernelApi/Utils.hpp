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
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/Mesh1D.hpp>
#include <MeshKernel/Mesh2D.hpp>
#include <MeshKernel/MeshEdgeCenters.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/Splines.hpp>

#include <MeshKernelApi/GeometryList.hpp>
#include <MeshKernelApi/GriddedSamples.hpp>
#include <MeshKernelApi/Mesh1D.hpp>
#include <MeshKernelApi/Mesh2D.hpp>

#include "MeshKernel/BilinearInterpolationOnGriddedSamples.hpp"
#include <MeshKernel/CurvilinearGrid/CurvilinearGridRectangular.hpp>

#include "MeshKernelApi/CurvilinearGrid.hpp"

#include <optional>
#include <span>
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

    /// @brief Computes the samples represented in gridded data in a vector of samples
    /// @param[in] griddedSamples The gridded data to convert
    /// @returns The converted vector of samples
    template <meshkernel::InterpolatableType T>
    static std::vector<meshkernel::Sample> ComputeGriddedDataSamples(const GriddedSamples& griddedSamples)
    {
        std::vector<meshkernel::Sample> result;
        meshkernel::Point origin{griddedSamples.x_origin, griddedSamples.y_origin};
        const auto numSamples = static_cast<size_t>(griddedSamples.num_x * griddedSamples.num_y);
        result.resize(numSamples);
        const T* valuePtr = static_cast<T*>(griddedSamples.values);
        if (griddedSamples.x_coordinates == nullptr || griddedSamples.y_coordinates == nullptr)
        {
            meshkernel::UInt index = 0;

            for (int j = 0; j < griddedSamples.num_x; ++j)
            {
                for (int i = griddedSamples.num_y - 1; i >= 0; --i)
                {
                    const auto griddedIndex = griddedSamples.num_x * i + j;
                    result[index].x = origin.x + j * griddedSamples.cell_size;
                    result[index].y = origin.y + i * griddedSamples.cell_size;
                    result[index].value = static_cast<double>(valuePtr[griddedIndex]);
                    index++;
                }
            }
            return result;
        }

        meshkernel::UInt index = 0;
        for (int j = 0; j < griddedSamples.num_x; ++j)
        {
            for (int i = griddedSamples.num_y - 1; i >= 0; --i)
            {
                const auto griddedIndex = griddedSamples.num_x * i + j;
                result[index].x = origin.x + griddedSamples.x_coordinates[griddedIndex];
                result[index].y = origin.y + griddedSamples.y_coordinates[griddedIndex];
                result[index].value = static_cast<double>(valuePtr[griddedIndex]);
                index++;
            }
        }
        return result;
    }

    /// @brief Converts the samples represented in gridded data in a vector of samples
    /// @param[in] griddedSamples The gridded data to convert
    /// @returns The converted vector of samples
    static std::vector<meshkernel::Sample> ConvertGriddedData(const GriddedSamples& griddedSamples)
    {
        std::vector<meshkernel::Sample> result;
        if (griddedSamples.num_x <= 0 || griddedSamples.num_y <= 0)
        {
            return result;
        }

        if (griddedSamples.value_type == static_cast<int>(meshkernel::InterpolationDataTypes::Short))
        {
            return ComputeGriddedDataSamples<short>(griddedSamples);
        }
        if (griddedSamples.value_type == static_cast<int>(meshkernel::InterpolationDataTypes::Float))
        {
            return ComputeGriddedDataSamples<float>(griddedSamples);
        }
        if (griddedSamples.value_type == static_cast<int>(meshkernel::InterpolationDataTypes::Double))
        {
            return ComputeGriddedDataSamples<double>(griddedSamples);
        }
        if (griddedSamples.value_type == static_cast<int>(meshkernel::InterpolationDataTypes::Int))
        {
            return ComputeGriddedDataSamples<int>(griddedSamples);
        }
        throw meshkernel::MeshKernelError("The value type for the gridded data samples is invalid.");
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

        const auto indices = FindIndices(splineCornerPoints,
                                         0,
                                         static_cast<meshkernel::UInt>(splineCornerPoints.size()),
                                         meshkernel::constants::missing::doubleValue);

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
    static void SetMesh2dApiDimensions(const meshkernel::Mesh& mesh2d, Mesh2D& mesh2dApi)
    {
        size_t num_face_nodes = 0;
        for (size_t f = 0; f < mesh2d.GetNumFaces(); f++)
        {
            num_face_nodes += mesh2d.m_facesNodes[f].size();
        }

        mesh2dApi.num_face_nodes = static_cast<int>(num_face_nodes);
        mesh2dApi.num_faces = static_cast<int>(mesh2d.GetNumFaces());
        mesh2dApi.num_nodes = static_cast<int>(mesh2d.GetNumNodes());
        mesh2dApi.num_valid_nodes = static_cast<int>(mesh2d.GetNumValidNodes());
        mesh2dApi.num_edges = static_cast<int>(mesh2d.GetNumEdges());
        mesh2dApi.num_valid_edges = static_cast<int>(mesh2d.GetNumValidEdges());
    }

    /// @brief Sets only the node and edge data for meshkernelapi::Mesh2D
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetMesh2dApiNodeEdgeData(meshkernel::Mesh2D& mesh2d, Mesh2D& mesh2dApi)
    {
        if (mesh2dApi.num_nodes != static_cast<int>(mesh2d.GetNumNodes()))
        {
            throw meshkernel::ConstraintError("The number of nodes in the mesh2d api structure does not equal the number of nodes in the grid, {} /= {}",
                                              mesh2dApi.num_nodes, mesh2d.GetNumNodes());
        }

        if (mesh2dApi.num_edges != static_cast<int>(mesh2d.GetNumEdges()))
        {
            throw meshkernel::ConstraintError("The number of edges in the mesh2d api structure does not equal the number of edges in the grid, {} /= {}",
                                              mesh2dApi.num_edges, mesh2d.GetNumEdges());
        }

        for (meshkernel::UInt n = 0; n < mesh2d.GetNumNodes(); ++n)
        {
            mesh2dApi.node_x[n] = mesh2d.Node(n).x;
            mesh2dApi.node_y[n] = mesh2d.Node(n).y;
        }

        for (meshkernel::UInt edgeIndex = 0; edgeIndex < mesh2d.GetNumEdges(); ++edgeIndex)
        {
            mesh2dApi.edge_nodes[edgeIndex * 2] = static_cast<int>(mesh2d.GetEdge(edgeIndex).first);
            mesh2dApi.edge_nodes[edgeIndex * 2 + 1] = static_cast<int>(mesh2d.GetEdge(edgeIndex).second);
        }
        SetMesh2dApiDimensions(mesh2d, mesh2dApi);
    }

    /// @brief Sets the meshkernelapi::Mesh2D data
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetMesh2dApiData(meshkernel::Mesh2D& mesh2d, Mesh2D& mesh2dApi)
    {
        std::vector<meshkernel::Point> edgeCentres = meshkernel::algo::ComputeEdgeCentres(mesh2d);

        for (meshkernel::UInt n = 0; n < mesh2d.GetNumNodes(); n++)
        {
            mesh2dApi.node_x[n] = mesh2d.Node(n).x;
            mesh2dApi.node_y[n] = mesh2d.Node(n).y;
        }

        for (meshkernel::UInt edgeIndex = 0; edgeIndex < mesh2d.GetNumEdges(); edgeIndex++)
        {
            mesh2dApi.edge_x[edgeIndex] = edgeCentres[edgeIndex].x;
            mesh2dApi.edge_y[edgeIndex] = edgeCentres[edgeIndex].y;
            mesh2dApi.edge_nodes[edgeIndex * 2] = static_cast<int>(mesh2d.GetEdge(edgeIndex).first);
            mesh2dApi.edge_nodes[edgeIndex * 2 + 1] = static_cast<int>(mesh2d.GetEdge(edgeIndex).second);

            const auto& firstEdgeFace = mesh2d.m_edgesFaces[edgeIndex][0];
            mesh2dApi.edge_faces[edgeIndex * 2] = firstEdgeFace == meshkernel::constants::missing::uintValue ? -1 : static_cast<int>(firstEdgeFace);
            const auto& secondEdgeFace = mesh2d.m_edgesFaces[edgeIndex][1];
            mesh2dApi.edge_faces[edgeIndex * 2 + 1] = secondEdgeFace == meshkernel::constants::missing::uintValue ? -1 : static_cast<int>(secondEdgeFace);
        }

        int faceIndex = 0;
        for (size_t f = 0; f < mesh2d.GetNumFaces(); f++)
        {
            mesh2dApi.face_x[f] = mesh2d.m_facesMassCenters[f].x;
            mesh2dApi.face_y[f] = mesh2d.m_facesMassCenters[f].y;
            mesh2dApi.nodes_per_face[f] = static_cast<int>(mesh2d.m_facesNodes[f].size());
            for (size_t n = 0; n < mesh2d.m_facesNodes[f].size(); ++n)
            {
                mesh2dApi.face_nodes[faceIndex] = static_cast<int>(mesh2d.m_facesNodes[f][n]);
                mesh2dApi.face_edges[faceIndex] = static_cast<int>(mesh2d.m_facesEdges[f][n]);
                faceIndex++;
            }
        }
        SetMesh2dApiDimensions(mesh2d, mesh2dApi);
    }

    /// @brief Sets the meshkernelapi::CurvilinearGrid data
    /// @param[in]  mesh2d    The meshkernel::Mesh2D instance
    /// @param[out] mesh2dApi The output meshkernelapi::Mesh2D instance
    static void SetCurvilinearGridApiData(const meshkernel::CurvilinearGrid& curvilinearGrid,
                                          CurvilinearGrid& curvilinearGridApi)
    {
        if (curvilinearGridApi.num_n != static_cast<int>(curvilinearGrid.NumN()))
        {
            throw meshkernel::ConstraintError("The number of rows in the api structure does not equal the number of rows in the grid, {} /= {}",
                                              curvilinearGridApi.num_n, curvilinearGrid.NumN());
        }

        if (curvilinearGridApi.num_m != static_cast<int>(curvilinearGrid.NumM()))
        {
            throw meshkernel::ConstraintError("The number of columns in the api structure does not equal the number of columns in the grid, {} /= {}",
                                              curvilinearGridApi.num_m, curvilinearGrid.NumM());
        }

        int count = 0;

        for (meshkernel::UInt n = 0; n < curvilinearGrid.NumN(); n++)
        {
            for (meshkernel::UInt m = 0; m < curvilinearGrid.NumM(); m++)
            {
                curvilinearGridApi.node_x[count] = curvilinearGrid.GetNode(n, m).x;
                curvilinearGridApi.node_y[count] = curvilinearGrid.GetNode(n, m).y;
                ++count;
            }
        }
    }

    /// @brief Sets dimensions members of meshkernelapi::Mesh1D instance
    /// @param[in]  mesh1d    The meshkernel::Mesh1D instance
    /// @param[out] mesh1dApi The output meshkernelapi::Mesh1D instance
    static void SetMesh1dApiDimension(const meshkernel::Mesh1D& mesh1d,
                                      Mesh1D& mesh1dApi)
    {
        mesh1dApi.num_nodes = static_cast<int>(mesh1d.GetNumNodes());
        mesh1dApi.num_valid_nodes = static_cast<int>(mesh1d.GetNumValidNodes());
        mesh1dApi.num_edges = static_cast<int>(mesh1d.GetNumEdges());
        mesh1dApi.num_valid_edges = static_cast<int>(mesh1d.GetNumValidEdges());
    }

    /// @brief Sets a meshkernelapi::Mesh1D data
    /// @param[in]  mesh1d           The input meshkernel::Mesh1D instance
    /// @param[out] mesh1dApi        The output meshkernelapi::Mesh1D instance
    static void SetMesh1dApiData(const meshkernel::Mesh1D& mesh1d,
                                 Mesh1D& mesh1dApi)
    {
        for (meshkernel::UInt n = 0; n < mesh1d.GetNumNodes(); n++)
        {
            mesh1dApi.node_x[n] = mesh1d.Node(n).x;
            mesh1dApi.node_y[n] = mesh1d.Node(n).y;
        }

        size_t edgeIndex = 0;
        for (meshkernel::UInt e = 0; e < mesh1d.GetNumEdges(); e++)
        {
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d.GetEdge(e).first);
            edgeIndex++;
            mesh1dApi.edge_nodes[edgeIndex] = static_cast<int>(mesh1d.GetEdge(e).second);
            edgeIndex++;
        }
        SetMesh1dApiDimension(mesh1d, mesh1dApi);
    }

    /// @brief Generate a rectangular curvilinear grid
    /// @param[in] makeGridParameters The parameters for creating a rectangular curvilinear grid are as follows
    /// @param[in] projection         The projection tu use
    /// @returns The generated curvilinear grid
    static std::unique_ptr<meshkernel::CurvilinearGrid> CreateRectangularCurvilinearGrid(const meshkernel::MakeGridParameters& makeGridParameters,
                                                                                         const meshkernel::Projection& projection)
    {
        meshkernel::CurvilinearGridRectangular grid(projection);
        return grid.Compute(makeGridParameters.num_columns,
                            makeGridParameters.num_rows,
                            makeGridParameters.origin_x,
                            makeGridParameters.origin_y,
                            makeGridParameters.angle,
                            makeGridParameters.block_size_x,
                            makeGridParameters.block_size_y);
    }

    /// @brief Generate a rectangular curvilinear grid from polygons
    /// @param[in] makeGridParameters The parameters for creating a rectangular curvilinear grid are as follows
    /// @param[in] geometryList       The polygon inside which generating the curvilinear grid
    /// @param[in] projection         The projection tu use
    /// @returns The generated curvilinear grid
    static std::unique_ptr<meshkernel::CurvilinearGrid> CreateRectangularCurvilinearGridFromPolygons(const meshkernel::MakeGridParameters& makeGridParameters,
                                                                                                     const GeometryList& geometryList,
                                                                                                     const meshkernel::Projection& projection)
    {
        const meshkernel::CurvilinearGridRectangular grid(projection);

        auto polygonNodes = ConvertGeometryListToPointVector(geometryList);

        const auto polygon = std::make_shared<meshkernel::Polygons>(polygonNodes, projection);

        return grid.Compute(makeGridParameters.angle,
                            makeGridParameters.block_size_x,
                            makeGridParameters.block_size_y,
                            polygon,
                            0);
    }

    /// @brief Generate a rectangular curvilinear grid based on extension
    /// @param[in] makeGridParameters The parameters for creating a rectangular curvilinear grid are as follows
    /// @param[in] projection         The projection tu use
    /// @returns The generated curvilinear grid
    static std::unique_ptr<meshkernel::CurvilinearGrid> CreateRectangularCurvilinearGridOnExtension(const meshkernel::MakeGridParameters& makeGridParameters,
                                                                                                    const meshkernel::Projection& projection)
    {
        const meshkernel::CurvilinearGridRectangular grid(projection);

        if (!meshkernel::IsEqual(makeGridParameters.angle, 0.0))
        {
            throw meshkernel::AlgorithmError("When generating an uniform grid on an defined extension, the grid angle must be equal to 0");
        }

        return grid.Compute(makeGridParameters.origin_x,
                            makeGridParameters.origin_y,
                            makeGridParameters.block_size_x,
                            makeGridParameters.block_size_y,
                            makeGridParameters.upper_right_x,
                            makeGridParameters.upper_right_y);
    }

    template <meshkernel::InterpolatableType T>
    static std::unique_ptr<meshkernel::MeshInterpolation> CreateBilinearInterpolator(const meshkernel::Mesh2D& mesh2d,
                                                                                     const GriddedSamples& griddedSamples)
    {
        meshkernel::Point origin{griddedSamples.x_origin, griddedSamples.y_origin};
        if (griddedSamples.x_coordinates == nullptr || griddedSamples.y_coordinates == nullptr)
        {
            return std::make_unique<meshkernel::BilinearInterpolationOnGriddedSamples<T>>(mesh2d,
                                                                                          griddedSamples.num_x,
                                                                                          griddedSamples.num_y,
                                                                                          origin,
                                                                                          griddedSamples.cell_size,
                                                                                          std::span<T const>{reinterpret_cast<T const* const>(griddedSamples.values),
                                                                                                             static_cast<size_t>(griddedSamples.num_x * griddedSamples.num_y)});
        }
        return std::make_unique<meshkernel::BilinearInterpolationOnGriddedSamples<T>>(mesh2d,
                                                                                      std::span<double const>{griddedSamples.x_coordinates,
                                                                                                              static_cast<size_t>(griddedSamples.num_x)},
                                                                                      std::span<double const>{griddedSamples.y_coordinates,
                                                                                                              static_cast<size_t>(griddedSamples.num_y)},
                                                                                      std::span<T const>{reinterpret_cast<T const* const>(griddedSamples.values),
                                                                                                         static_cast<size_t>(griddedSamples.num_x * griddedSamples.num_y)});
    }

    static std::unique_ptr<meshkernel::MeshInterpolation> CreateBilinearInterpolatorBasedOnType(const GriddedSamples& griddedSamples,
                                                                                                const meshkernel::Mesh2D& mesh2d)
    {

        switch (static_cast<meshkernel::InterpolationDataTypes>(griddedSamples.value_type))
        {
        case meshkernel::InterpolationDataTypes::Short:
            return CreateBilinearInterpolator<short>(mesh2d, griddedSamples);
        case meshkernel::InterpolationDataTypes::Float:
            return CreateBilinearInterpolator<float>(mesh2d, griddedSamples);
        case meshkernel::InterpolationDataTypes::Int:
            return CreateBilinearInterpolator<int>(mesh2d, griddedSamples);
        case meshkernel::InterpolationDataTypes::Double:
            return CreateBilinearInterpolator<double>(mesh2d, griddedSamples);
        default:
            throw meshkernel::MeshKernelError("Invalid value_type for GriddedSamples");
        }
    }

    static void FillFacePolygons(const meshkernel::Mesh2D& mesh2d,
                                 const std::vector<bool>& facesInPolygon,
                                 const GeometryList& facePolygons)
    {
        meshkernel::UInt count = 0;
        for (meshkernel::UInt f = 0; f < mesh2d.GetNumFaces(); ++f)
        {
            if (!facesInPolygon[f])
            {
                continue;
            }

            const auto& faceNodes = mesh2d.m_facesNodes[f];
            if (count != 0)
            {
                facePolygons.coordinates_x[count] = meshkernel::constants::missing::doubleValue;
                facePolygons.coordinates_y[count] = meshkernel::constants::missing::doubleValue;
                count++;
            }

            for (meshkernel::UInt n = 0u; n < faceNodes.size(); ++n)
            {
                const auto& currentNode = mesh2d.Node(faceNodes[n]);
                facePolygons.coordinates_x[count] = currentNode.x;
                facePolygons.coordinates_y[count] = currentNode.y;
                count++;
            }
            const auto& currentNode = mesh2d.Node(faceNodes[0]);
            facePolygons.coordinates_x[count] = currentNode.x;
            facePolygons.coordinates_y[count] = currentNode.y;
            count++;
        }
    }

    static std::optional<double> MinValidElement(const std::vector<double>& values)
    {
        auto filtered = values | std::views::filter([](const double& v)
                                                    { return v != meshkernel::constants::missing::doubleValue; });
        auto begin = std::ranges::begin(filtered);
        auto end = std::ranges::end(filtered);

        if (begin == end)
        {
            return std::nullopt; // No valid elements
        }

        return *std::ranges::min_element(filtered);
    }

} // namespace meshkernelapi
