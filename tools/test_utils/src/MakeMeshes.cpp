#include <array>
#include <memory>
#include <stdexcept>

#include <netcdf.h>

#include <TestUtils/MakeMeshes.hpp>

std::tuple<size_t,
           size_t,
           std::vector<double>,
           std::vector<double>,
           std::vector<int>,
           std::vector<int>,
           std::vector<int>>
ReadLegacyMeshFile(std::filesystem::path const& file_path)
{
    {
        std::error_code error_code;
        if (!std::filesystem::exists(file_path, error_code))
        {
            throw std::filesystem::filesystem_error("File does not exist", error_code);
        }
    }

    int ncidp = 0;
    int err = nc_open(file_path.string().c_str(), NC_NOWRITE, &ncidp);
    if (err != 0)
    {
        throw "ReadLegacyMesh2DFromFile: Could not load netcdf file.";
    }

    std::string meshNodesName{"nNetNode"};
    int dimid = 0;
    err = nc_inq_dimid(ncidp, meshNodesName.c_str(), &dimid);
    if (err != 0)
    {
        throw "ReadLegacyMesh2DFromFile: Could not find the ID of a dimension of 'nNetNode'.";
    }

    std::size_t num_nodes;
    std::array<char, NC_MAX_NAME> read_name;
    read_name.fill('\0');
    err = nc_inq_dim(ncidp, dimid, read_name.data(), &num_nodes);
    if (err != 0)
    {
        throw "ReadLegacyMesh2DFromFile: Could not find the length of dimension of 'nNetNode'.";
    }

    std::string meshEdgesName{"nNetLink"};
    err = nc_inq_dimid(ncidp, meshEdgesName.c_str(), &dimid);
    if (err != 0)
    {
        throw "ReadLegacyMesh2DFromFile: Could not find the ID of a dimension of 'nNetLink'.";
    }

    std::size_t num_edges;
    nc_inq_dim(ncidp, dimid, read_name.data(), &num_edges);

    std::vector<double> node_x(num_nodes);
    std::vector<double> node_y(num_nodes);
    std::vector<int> edge_nodes(num_edges * 2);
    std::vector<int> edge_type(num_edges);

    std::string meshNodeXName{"NetNode_x"};
    int varid = 0;
    nc_inq_varid(ncidp, meshNodeXName.c_str(), &varid);
    nc_get_var_double(ncidp, varid, node_x.data());

    std::string meshNodeYName{"NetNode_y"};
    nc_inq_varid(ncidp, meshNodeYName.c_str(), &varid);
    nc_get_var_double(ncidp, varid, node_y.data());

    std::string linkName{"NetLink"};
    nc_inq_varid(ncidp, linkName.c_str(), &varid);
    nc_get_var_int(ncidp, varid, edge_nodes.data());

    std::string edgeTypeName{"NetLinkType"};
    nc_inq_varid(ncidp, edgeTypeName.c_str(), &varid);
    nc_get_var_int(ncidp, varid, edge_type.data());

    nc_close(ncidp);

    // Transform into 0 based indexing
    for (size_t i = 0; i < num_edges * 2; i++)
    {
        edge_nodes[i] -= 1;
    }

    std::vector<int> node_type(num_nodes);
    size_t index = 0;
    for (size_t i = 0; i < num_edges; i++)
    {
        const auto type = edge_type[i];
        const auto firstNode = edge_nodes[index];
        index++;
        const auto secondNode = edge_nodes[index];
        index++;
        node_type[firstNode] = type;
        node_type[secondNode] = type;
    }

    return {num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type};
}

std::tuple<std::vector<meshkernel::Point>,
           std::vector<meshkernel::Edge>>
ComputeEdgesAndNodes(
    std::filesystem::path const& file_path,
    meshkernel::Mesh::Type meshType)
{
    const auto [num_nodes, num_edges, node_x, node_y, node_type, edge_nodes, edge_type] =
        ReadLegacyMeshFile(file_path);
    std::vector<meshkernel::Edge> edges;
    std::vector<meshkernel::Point> nodes;
    std::vector<meshkernel::UInt> nodeMapping;
    edges.reserve(num_edges);
    nodes.reserve(num_nodes);
    nodeMapping.resize(num_nodes);

    int nodeType = 0;
    if (meshType == meshkernel::Mesh::Type::Mesh1D)
    {
        nodeType = 1;
    }
    if (meshType == meshkernel::Mesh::Type::Mesh2D)
    {
        nodeType = 2;
    }

    for (meshkernel::UInt i = 0; i < num_nodes; i++)
    {
        // If the node is not part of a 2 mesh, do not add it in nodes
        if (node_type[i] == nodeType || node_type[i] == 0)
        {
            nodes.emplace_back(node_x[i], node_y[i]);
            nodeMapping[i] = static_cast<meshkernel::UInt>(nodes.size()) - 1;
        }
    }

    auto index = 0;
    for (meshkernel::UInt i = 0; i < num_edges; i++)
    {

        auto const firstNode = edge_nodes[index];
        auto const secondNode = edge_nodes[index + 1];
        if ((node_type[firstNode] == nodeType || node_type[firstNode] == 0) &&
            (node_type[secondNode] == nodeType || node_type[firstNode] == 0))
        {
            edges.emplace_back(nodeMapping[firstNode], nodeMapping[secondNode]);
        }
        index = index + 2;
    }

    return {nodes, edges};
}

std::shared_ptr<meshkernel::Mesh2D> ReadLegacyMesh2DFromFile(
    std::filesystem::path const& file_path,
    meshkernel::Projection projection)
{
    const auto [nodes, edges] = ComputeEdgesAndNodes(file_path, meshkernel::Mesh::Type::Mesh2D);
    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, projection);
}

std::shared_ptr<meshkernel::Mesh1D> ReadLegacyMesh1DFromFile(
    std::filesystem::path const& file_path,
    meshkernel::Projection projection)
{
    const auto [nodes, edges] = ComputeEdgesAndNodes(file_path, meshkernel::Mesh::Type::Mesh1D);
    return std::make_shared<meshkernel::Mesh1D>(edges, nodes, projection);
}

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(
    meshkernel::UInt n,
    meshkernel::UInt m,
    double dim_x,
    double dim_y,
    meshkernel::Projection projection,
    meshkernel::Point const& origin)
{
    std::vector<std::vector<meshkernel::UInt>> node_indices(n, std::vector<meshkernel::UInt>(m));
    std::vector<meshkernel::Point> nodes(n * m);

    {
        meshkernel::UInt index = 0;
        double const delta_x = dim_x / static_cast<double>(n - 1);
        double const delta_y = dim_y / static_cast<double>(m - 1);
        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            for (meshkernel::UInt j = 0; j < m; ++j)
            {
                node_indices[i][j] = i * m + j;
                nodes[index] = {origin.x + i * delta_x, origin.y + j * delta_y};
                index++;
            }
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);

    {
        meshkernel::UInt index = 0;

        for (meshkernel::UInt i = 0; i < n - 1; ++i)
        {
            for (meshkernel::UInt j = 0; j < m; ++j)
            {
                edges[index] = {node_indices[i][j], node_indices[i + 1][j]};
                index++;
            }
        }

        for (meshkernel::UInt i = 0; i < n; ++i)
        {
            for (meshkernel::UInt j = 0; j < m - 1; ++j)
            {
                edges[index] = {node_indices[i][j + 1], node_indices[i][j]};
                index++;
            }
        }
    }

    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, projection);
}

std::shared_ptr<meshkernel::Mesh2D> MakeRectangularMeshForTesting(
    meshkernel::UInt n,
    meshkernel::UInt m,
    double delta,
    meshkernel::Projection projection,
    meshkernel::Point const& origin)
{
    double const dim_x = delta * static_cast<double>(n - 1);
    double const dim_y = delta * static_cast<double>(m - 1);
    return MakeRectangularMeshForTesting(
        n,
        m,
        dim_x,
        dim_y,
        projection,
        origin);
}

std::tuple<meshkernel::UInt,
           meshkernel::UInt,
           std::vector<double>,
           std::vector<double>,
           std::vector<int>>
MakeRectangularMeshForApiTesting(
    meshkernel::UInt numRows,
    meshkernel::UInt numColumns,
    double delta)
{

    const auto numY = numRows + static_cast<meshkernel::UInt>(1);
    const auto numX = numColumns + static_cast<meshkernel::UInt>(1);

    std::vector<std::vector<meshkernel::UInt>> indicesValues(numX, std::vector<meshkernel::UInt>(numY));

    std::vector<double> nodeX(numX * numY);
    std::vector<double> nodeY(numX * numY);

    meshkernel::UInt nodeIndex = 0;
    for (meshkernel::UInt i = 0; i < numX; ++i)
    {
        for (meshkernel::UInt j = 0; j < numY; ++j)
        {

            nodeX[nodeIndex] = i * delta;
            nodeY[nodeIndex] = j * delta;
            indicesValues[i][j] = i * numY + j;
            nodeIndex++;
        }
    }

    std::vector<int> edgeNodes(((numX - 1) * numY + (numY - 1) * numX) * 2);
    meshkernel::UInt edgeIndex = 0;
    for (meshkernel::UInt i = 0; i < numX - 1; ++i)
    {
        for (meshkernel::UInt j = 0; j < numY; ++j)
        {
            edgeNodes[edgeIndex] = static_cast<int>(indicesValues[i][j]);
            edgeIndex++;
            edgeNodes[edgeIndex] = static_cast<int>(indicesValues[i + 1][j]);
            edgeIndex++;
        }
    }

    for (meshkernel::UInt i = 0; i < numX; ++i)
    {
        for (meshkernel::UInt j = 0; j < numY - 1; ++j)
        {
            edgeNodes[edgeIndex] = static_cast<int>(indicesValues[i][j + 1]);
            edgeIndex++;
            edgeNodes[edgeIndex] = static_cast<int>(indicesValues[i][j]);
            edgeIndex++;
        }
    }

    const meshkernel::UInt numNodes = nodeIndex;
    const meshkernel::UInt numEdges = static_cast<meshkernel::UInt>(edgeIndex * 0.5);

    return {numNodes, numEdges, nodeX, nodeY, edgeNodes};
}

std::shared_ptr<meshkernel::Mesh2D> MakeSmallSizeTriangularMeshForTestingAsNcFile()
{
    // Prepare
    std::vector<meshkernel::Point> nodes;

    nodes.push_back({322.252624511719, 454.880187988281});
    nodes.push_back({227.002044677734, 360.379241943359});
    nodes.push_back({259.252227783203, 241.878051757813});
    nodes.push_back({428.003295898438, 210.377746582031});
    nodes.push_back({536.003967285156, 310.878753662109});
    nodes.push_back({503.753784179688, 432.379974365234});
    nodes.push_back({350.752807617188, 458.630249023438});
    nodes.push_back({343.15053976393, 406.232256102912});
    nodes.push_back({310.300984548069, 319.41005739802});
    nodes.push_back({423.569603308318, 326.17986967523});

    std::vector<meshkernel::Edge> edges;
    edges.push_back({2, 8});
    edges.push_back({1, 8});
    edges.push_back({1, 2});
    edges.push_back({2, 3});
    edges.push_back({3, 8});
    edges.push_back({1, 7});
    edges.push_back({0, 7});
    edges.push_back({0, 1});
    edges.push_back({7, 8});
    edges.push_back({6, 7});
    edges.push_back({0, 6});
    edges.push_back({8, 9});
    edges.push_back({7, 9});
    edges.push_back({3, 4});
    edges.push_back({4, 9});
    edges.push_back({3, 9});
    edges.push_back({5, 7});
    edges.push_back({5, 6});
    edges.push_back({5, 9});
    edges.push_back({4, 5});

    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);
}

std::shared_ptr<meshkernel::Mesh2D> MakeCurvilinearGridForTesting()
{
    std::vector<double> xCoordinates{777.2400642395231,
                                     776.8947176796199,
                                     776.5500495969297,
                                     776.2062753740686,
                                     775.8639152055595,
                                     789.0362557653892,
                                     788.6869233110746,
                                     788.3382299013459,
                                     787.9903970612199,
                                     787.6439255320029,
                                     800.8316506843166,
                                     800.4781249366497,
                                     800.125196925116,
                                     799.7730958077465,
                                     799.4223012497744,
                                     783.1381600024562,
                                     782.7908204953472,
                                     782.4441397491378,
                                     782.0983362176443,
                                     781.7539203687811,
                                     794.933953224853,
                                     794.5825241238622,
                                     794.231713413231,
                                     793.8817464344831,
                                     793.5331133908886,
                                     777.0673909595715,
                                     776.7223836382748,
                                     776.3781624854992,
                                     776.0350952898141,
                                     788.8615895382319,
                                     788.5125766062102,
                                     788.164313481283,
                                     787.8171612966114,
                                     800.6548878104832,
                                     800.3016609308829,
                                     799.9491463664313,
                                     799.5976985287605,
                                     782.9644902489017,
                                     782.6174801222426,
                                     782.2712379833911,
                                     781.9261282932127,
                                     794.7582386743576,
                                     794.4071187685466,
                                     794.0567299238571,
                                     793.7074299126859};

    std::vector<double> yCoordinates{1145.6125094106028,
                                     1147.65972166567,
                                     1149.7057403812394,
                                     1151.7507901691963,
                                     1153.7936194659985,
                                     1147.830596354202,
                                     1149.8579007347682,
                                     1151.884213062344,
                                     1153.9097092742677,
                                     1155.9332418178237,
                                     1150.0482786995756,
                                     1152.0565093840655,
                                     1154.063958611566,
                                     1156.0707463059343,
                                     1158.0758372662724,
                                     1146.7215528824024,
                                     1148.758811200219,
                                     1150.7949767217917,
                                     1152.830249721732,
                                     1154.8634306419112,
                                     1148.939437526889,
                                     1150.957205059417,
                                     1152.974085836955,
                                     1154.990227790101,
                                     1157.004539542048,
                                     1146.6361155381364,
                                     1148.6827310234548,
                                     1150.7282652752178,
                                     1152.7722048175974,
                                     1148.8442485444853,
                                     1150.8710568985562,
                                     1152.896961168306,
                                     1154.9214755460457,
                                     1151.0523940418207,
                                     1153.0602339978159,
                                     1155.06735245875,
                                     1157.0732917861033,
                                     1147.7401820413108,
                                     1149.7768939610055,
                                     1151.812613221762,
                                     1153.8468401818216,
                                     1149.948321293153,
                                     1151.965645448186,
                                     1153.982156813528,
                                     1155.9973836660745};

    std::vector<meshkernel::Edge> edges{{1, 16},
                                        {2, 17},
                                        {3, 18},
                                        {4, 19},
                                        {5, 20},
                                        {6, 21},
                                        {7, 22},
                                        {8, 23},
                                        {9, 24},
                                        {10, 25},
                                        {1, 26},
                                        {2, 27},
                                        {3, 28},
                                        {4, 29},
                                        {6, 30},
                                        {7, 31},
                                        {8, 32},
                                        {9, 33},
                                        {11, 34},
                                        {12, 35},
                                        {13, 36},
                                        {14, 37},
                                        {16, 38},
                                        {30, 38},
                                        {17, 38},
                                        {26, 38},
                                        {17, 39},
                                        {31, 39},
                                        {18, 39},
                                        {27, 39},
                                        {18, 40},
                                        {32, 40},
                                        {19, 40},
                                        {28, 40},
                                        {19, 41},
                                        {33, 41},
                                        {20, 41},
                                        {29, 41},
                                        {21, 42},
                                        {34, 42},
                                        {22, 42},
                                        {30, 42},
                                        {22, 43},
                                        {35, 43},
                                        {23, 43},
                                        {31, 43},
                                        {23, 44},
                                        {36, 44},
                                        {24, 44},
                                        {32, 44},
                                        {24, 45},
                                        {37, 45},
                                        {25, 45},
                                        {33, 45},
                                        {6, 16},
                                        {7, 17},
                                        {8, 18},
                                        {9, 19},
                                        {10, 20},
                                        {11, 21},
                                        {12, 22},
                                        {13, 23},
                                        {14, 24},
                                        {15, 25},
                                        {2, 26},
                                        {3, 27},
                                        {4, 28},
                                        {5, 29},
                                        {7, 30},
                                        {8, 31},
                                        {9, 32},
                                        {10, 33},
                                        {12, 34},
                                        {13, 35},
                                        {14, 36},
                                        {15, 37}};

    std::vector<meshkernel::Point> nodes(xCoordinates.size());

    for (size_t i = 0; i < nodes.size(); i++)
    {
        nodes[i].x = xCoordinates[i];
        nodes[i].y = yCoordinates[i];
    }

    for (size_t i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }
    return std::make_shared<meshkernel::Mesh2D>(edges, nodes, meshkernel::Projection::cartesian);
}

/// @brief Make face nodes
std::tuple<std::vector<double>,
           std::vector<double>,
           std::vector<int>,
           std::vector<int>,
           std::vector<int>>
MakeMeshWithFaceNodes()
{
    // Set-up new mesh
    std::vector<double> nodes_x{
        0,
        10,
        20,
        30,
        0,
        10,
        20,
        30,
        0,
        10,
        20,
        30,
        0,
        10,
        20,
        30};

    std::vector<double> nodes_y{
        0,
        0,
        0,
        0,
        10,
        10,
        10,
        10,
        20,
        20,
        20,
        20,
        30,
        30,
        30,
        30};

    std::vector<int> edges{
        0,
        1,

        1,
        5,

        5,
        4,

        4,
        0,

        1,
        2,

        2,
        6,

        6,
        5,

        2,
        3,

        3,
        7,

        7,
        6,

        4,
        8,

        5,
        9,

        9,
        8,

        6,
        10,

        10,
        9,

        7,
        11,

        11,
        10,

        8,
        12,

        9,
        13,

        13,
        12,

        10,
        14,

        14,
        13,

        11,
        15,

        15,
        14};

    std::vector<int> faceNodes{
        0,
        1,
        5,
        4,

        1,
        2,
        6,
        5,

        2,
        3,
        7,
        6,

        4,
        5,
        9,
        8,

        5,
        6,
        10,
        9,

        6,
        7,
        11,
        10,

        8,
        9,
        13,
        12,

        9,
        10,
        14,
        13,

        10,
        11,
        15,
        14};

    std::vector<int> num_face_nodes{4, 4, 4, 4, 4, 4, 4, 4, 4};

    return {nodes_x, nodes_y, edges, faceNodes, num_face_nodes};
}
