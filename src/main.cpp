#include "Orthogonalization.hpp"
#include "Mesh.hpp"
#include "../thirdParty/netcdf/include/netcdf.h"
#include <Windows.h>
#include <chrono>
#include <string>
#include <iostream>

int main()
{
    const std::size_t ENV_BUF_SIZE = 2000;
    char buf[ENV_BUF_SIZE];
    std::size_t bufsize = ENV_BUF_SIZE;
    int err = getenv_s(&bufsize, buf, bufsize, "PATH");
    if (err)
    {
        exit(EXIT_FAILURE);
    }
    std::string env_path = buf;
    err = _putenv_s("PATH", env_path.c_str());
    if (err) {
        std::cerr << "`_putenv_s` failed, returned " << err << std::endl;
        exit(EXIT_FAILURE);
    }

    auto netcdf = LoadLibrary("netcdf.dll");
    if (!netcdf)
    {
        exit(EXIT_FAILURE);
    }

    //get the mesh from file
    std::string c_path{ "D:\\LUCA\\ENGINES\\GridGeom++\\lastCheckout\\tests\\TestOrthogonalizationLargeTriangularGrid_net.nc" };
    typedef  int(__stdcall* nc_open_dll)(const char *path, int mode, int *ncidp);
    auto nc_open = (nc_open_dll)GetProcAddress(netcdf, "nc_open");

    typedef  int(__stdcall* nc_inq_dimid_dll)(int ncid, const char *name, int *idp);
    auto nc_inq_dimid = (nc_inq_dimid_dll)GetProcAddress(netcdf, "nc_inq_dimid");

    typedef  int(__stdcall* nc_inq_dim_dll)(int ncid, int dimid, char *name, size_t *lenp);
    auto nc_inq_dim = (nc_inq_dim_dll)GetProcAddress(netcdf, "nc_inq_dim");

    typedef  int(__stdcall* nc_inq_varid_dll)(int ncid, const char *name, int *varidp);
    auto nc_inq_varid = (nc_inq_varid_dll)GetProcAddress(netcdf, "nc_inq_varid");

    typedef  int(__stdcall* nc_get_var_double_dll)(int ncid, int varid, double* ip);
    auto nc_get_var_double = (nc_get_var_double_dll)GetProcAddress(netcdf, "nc_get_var_double");

    typedef  int(__stdcall * nc_get_var_int_dll)(int ncid, int varid, int* ip);
    auto nc_get_var_int = (nc_get_var_int_dll)GetProcAddress(netcdf, "nc_get_var_int");

    int ncidp = 0;
    err = nc_open(c_path.c_str(), NC_NOWRITE, &ncidp);
 

    std::string mesh2dNodes{ "nNetNode" };
    int dimid = 0;
    err = nc_inq_dimid(ncidp, mesh2dNodes.c_str(), &dimid);

    size_t num_nodes;
    auto read_name = new char[NC_MAX_NAME];
    err = nc_inq_dim(ncidp, dimid, read_name, &num_nodes);


    std::string mesh2dEdges{ "nNetLink" };
    int dimidedges = 0;
    err = nc_inq_dimid(ncidp, mesh2dEdges.c_str(), &dimid);

    size_t num_edges;
    err = nc_inq_dim(ncidp, dimid, read_name, &num_edges);

    std::vector<double> nodeX(num_nodes, 0.0);
    std::vector<double> nodeY(num_nodes, 0.0);
    std::string mesh2dNodeX{ "NetNode_x" };
    int varid = 0;
    err = nc_inq_varid(ncidp, mesh2dNodeX.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeX[0]);

    std::string mesh2dNodeY{ "NetNode_y" };
    err = nc_inq_varid(ncidp, mesh2dNodeY.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeY[0]);


    std::string mesh2dEdgeNodes{ "NetLink" };
    err = nc_inq_varid(ncidp, mesh2dEdgeNodes.c_str(), &varid);

    std::vector<int> edge_nodes(num_edges * 2, 0.0);
    err = nc_get_var_int(ncidp, varid, &edge_nodes[0]);


    std::vector<GridGeom::Edge> edges(num_edges);
    std::vector<GridGeom::Point> nodes(num_nodes);

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

    // now build node-edge mapping
    
    std::cout << "start computing " << std::endl;
    auto start(std::chrono::steady_clock::now());
    GridGeom::Mesh mesh;
    mesh.setMesh(edges, nodes, GridGeom::Projections::cartesian);
    GridGeom::Orthogonalization orthogonalization;
    orthogonalization.initialize(mesh);
    orthogonalization.iterate(mesh);
    auto end(std::chrono::steady_clock::now());
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "Elapsed time " << duration << " s " << std::endl;
    std::cout << "Test finished " << std::endl;
	
	return 0;
}
