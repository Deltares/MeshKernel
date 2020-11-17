#pragma once

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif

#include <stdexcept>
#include <MeshKernel/Mesh.hpp>
#include <TestUtils/MakeMeshes.hpp>
#include "../../../extern/netcdf/netCDF 4.6.1/include/netcdf.h"

std::shared_ptr<meshkernel::Mesh> ReadLegacyMeshFromFile(std::string filePath, meshkernel::Projections projection)
{
    auto mesh = std::make_shared<meshkernel::Mesh>();

    auto netcdf = LoadLibrary("netcdf.dll");

    if (!netcdf)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not load 'netcdf.dll'.");
    }

    typedef int(__stdcall * nc_open_dll)(const char* path, int mode, int* ncidp);
    auto nc_open = (nc_open_dll)GetProcAddress(netcdf, "nc_open");

    typedef int(__stdcall * nc_inq_dimid_dll)(int ncid, const char* name, int* idp);
    auto nc_inq_dimid = (nc_inq_dimid_dll)GetProcAddress(netcdf, "nc_inq_dimid");

    typedef int(__stdcall * nc_inq_dim_dll)(int ncid, int dimid, char* name, std::size_t* lenp);
    auto nc_inq_dim = (nc_inq_dim_dll)GetProcAddress(netcdf, "nc_inq_dim");

    typedef int(__stdcall * nc_inq_varid_dll)(int ncid, const char* name, int* varidp);
    auto nc_inq_varid = (nc_inq_varid_dll)GetProcAddress(netcdf, "nc_inq_varid");

    typedef int(__stdcall * nc_get_var_double_dll)(int ncid, int varid, double* ip);
    auto nc_get_var_double = (nc_get_var_double_dll)GetProcAddress(netcdf, "nc_get_var_double");

    typedef int(__stdcall * nc_get_var_int_dll)(int ncid, int varid, int* ip);
    auto nc_get_var_int = (nc_get_var_int_dll)GetProcAddress(netcdf, "nc_get_var_int");

    int ncidp = 0;
    int err = nc_open(filePath.c_str(), NC_NOWRITE, &ncidp);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not load netcdf file.");
    }

    std::string mesh2dNodes{"nNetNode"};
    int dimid = 0;
    err = nc_inq_dimid(ncidp, mesh2dNodes.c_str(), &dimid);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not find the ID of a dimension of 'nNetNode'.");
    }

    std::size_t num_nodes;
    auto read_name = new char[NC_MAX_NAME];
    err = nc_inq_dim(ncidp, dimid, read_name, &num_nodes);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not gind the length of dimension of 'nNetNode'.");
    }

    std::string mesh2dEdges{"nNetLink"};
    err = nc_inq_dimid(ncidp, mesh2dEdges.c_str(), &dimid);
    if (err != 0)
    {
        throw std::invalid_argument("ReadLegacyMeshFromFile: Could not find the ID of a dimension of 'nNetLink'.");
    }

    std::size_t num_edges;
    err = nc_inq_dim(ncidp, dimid, read_name, &num_edges);
    delete[] read_name;

    std::vector<double> nodeX(num_nodes, 0.0);
    std::vector<double> nodeY(num_nodes, 0.0);
    std::vector<double> nodeZ(num_nodes, 0.0);
    std::string mesh2dNodeX{"NetNode_x"};
    int varid = 0;
    err = nc_inq_varid(ncidp, mesh2dNodeX.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeX[0]);

    std::string mesh2dNodeY{"NetNode_y"};
    err = nc_inq_varid(ncidp, mesh2dNodeY.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeY[0]);

    std::string mesh2dEdgeNodes{"NetLink"};
    err = nc_inq_varid(ncidp, mesh2dEdgeNodes.c_str(), &varid);

    std::vector<int> edge_nodes(num_edges * 2, 0);
    err = nc_get_var_int(ncidp, varid, &edge_nodes[0]);

    std::vector<meshkernel::Edge> edges(num_edges);
    std::vector<meshkernel::Point> nodes(num_nodes);

    for (int i = 0; i < nodeX.size(); i++)
    {
        nodes[i].x = nodeX[i];
        nodes[i].y = nodeY[i];
    }

    int index = 0;
    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first = edge_nodes[index] - 1;
        index++;

        edges[i].second = edge_nodes[index] - 1;
        index++;
    }

    mesh->Set(edges, nodes, projection);
    return mesh;
}

std::shared_ptr<meshkernel::Mesh> MakeSmallSizeTriangularMeshForTestingAsNcFile()
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

    auto mesh = std::make_shared<meshkernel::Mesh>();
    mesh->Set(edges, nodes, meshkernel::Projections::cartesian);

    return mesh;
}

std::shared_ptr<meshkernel::Mesh> MakeRectangularMeshForTesting(int n, int m, double delta, meshkernel::Projections projection, meshkernel::Point origin)
{
    std::vector<std::vector<int>> indexesValues(n, std::vector<int>(m));
    std::vector<meshkernel::Point> nodes(n * m);
    std::size_t nodeIndex = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            indexesValues[i][j] = i * m + j;
            nodes[nodeIndex] = {origin.x + i * delta, origin.y + j * delta};
            nodeIndex++;
        }
    }

    std::vector<meshkernel::Edge> edges((n - 1) * m + (m - 1) * n);
    std::size_t edgeIndex = 0;

    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            edges[edgeIndex] = {indexesValues[i][j], indexesValues[i + 1][j]};
            edgeIndex++;
        }
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m - 1; ++j)
        {
            edges[edgeIndex] = {indexesValues[i][j + 1], indexesValues[i][j]};
            edgeIndex++;
        }
    }

    auto mesh = std::make_shared<meshkernel::Mesh>();
    mesh->Set(edges, nodes, projection);

    return mesh;
}

std::shared_ptr<meshkernel::Mesh> MakeCurvilinearGridForTesting()
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

    for (int i = 0; i < nodes.size(); i++)
    {
        nodes[i].x = xCoordinates[i];
        nodes[i].y = yCoordinates[i];
    }

    for (int i = 0; i < edges.size(); i++)
    {
        edges[i].first -= 1;
        edges[i].second -= 1;
    }

    auto mesh = std::make_shared<meshkernel::Mesh>();
    mesh->Set(edges, nodes, meshkernel::Projections::cartesian);

    return mesh;
}
