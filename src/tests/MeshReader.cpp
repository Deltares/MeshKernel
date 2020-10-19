#pragma once

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#endif
#include <memory>
#include <iostream>

#include "../Mesh.hpp"
#include "../../thirdParty/third_party_open/netcdf/netCDF 4.6.1/include/netcdf.h"

static std::shared_ptr<meshkernel::Mesh> ReadLegacyMeshFromFile(std::string filePath)
{
    std::shared_ptr<meshkernel::Mesh> mesh = std::make_shared<meshkernel::Mesh>();

    auto netcdf = LoadLibrary("netcdf.dll");

    if (!netcdf)
    {
        return mesh;
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
        return mesh;
    }

    std::string mesh2dNodes{"nNetNode"};
    int dimid = 0;
    err = nc_inq_dimid(ncidp, mesh2dNodes.c_str(), &dimid);
    if (err != 0)
    {
        return mesh;
    }

    std::size_t num_nodes;
    auto read_name = new char[NC_MAX_NAME];
    err = nc_inq_dim(ncidp, dimid, read_name, &num_nodes);
    if (err != 0)
    {
        return mesh;
    }

    std::string mesh2dEdges{"nNetLink"};
    err = nc_inq_dimid(ncidp, mesh2dEdges.c_str(), &dimid);
    if (err != 0)
    {
        return mesh;
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

    mesh->Set(edges, nodes, meshkernel::Projections::cartesian);

    return mesh;
}