#include "Orthogonalization.hpp"
#include "Mesh.hpp"
#include "Entities.hpp"
#include "Gridgeom.hpp"
#include "../thirdParty/netcdf/include/netcdf.h"


#if defined(_WIN32)
#include <Windows.h>
#if defined (GRIDGEOM_API)
#undef GRIDGEOM_API
#define GRIDGEOM_API __declspec(dllimport)
#endif
#endif

#include <chrono>
#include <string>
#include <iostream>

using namespace GridGeomApi;

int main()
{
    std::cout << "load netcdf " << std::endl;
    auto netcdf = LoadLibrary("netcdf.dll");
    if (!netcdf)
    {
        exit(EXIT_FAILURE);
    }
    std::cout << "netcdf loaded " << netcdf << std::endl;

    //get the mesh from file
    std::string c_path{ "D:\\LUCA\\ENGINES\\GridGeom++\\lastCheckout\\tests\\TestOrthogonalizationLargeTriangularGrid_net.nc" };
    typedef  int(__stdcall* nc_open_dll)(const char *path, int mode, int *ncidp);
    auto nc_open = (nc_open_dll)GetProcAddress(netcdf, "nc_open");

    typedef  int(__stdcall* nc_inq_dimid_dll)(int ncid, const char *name, int *idp);
    auto nc_inq_dimid = (nc_inq_dimid_dll)GetProcAddress(netcdf, "nc_inq_dimid");

    typedef  int(__stdcall* nc_inq_dim_dll)(int ncid, int dimid, char *name, std::size_t *lenp);
    auto nc_inq_dim = (nc_inq_dim_dll)GetProcAddress(netcdf, "nc_inq_dim");

    typedef  int(__stdcall* nc_inq_varid_dll)(int ncid, const char *name, int *varidp);
    auto nc_inq_varid = (nc_inq_varid_dll)GetProcAddress(netcdf, "nc_inq_varid");

    typedef  int(__stdcall* nc_get_var_double_dll)(int ncid, int varid, double* ip);
    auto nc_get_var_double = (nc_get_var_double_dll)GetProcAddress(netcdf, "nc_get_var_double");

    typedef  int(__stdcall * nc_get_var_int_dll)(int ncid, int varid, int* ip);
    auto nc_get_var_int = (nc_get_var_int_dll)GetProcAddress(netcdf, "nc_get_var_int");

    int ncidp = 0;
    int err = nc_open(c_path.c_str(), NC_NOWRITE, &ncidp);

    std::string mesh2dNodes{ "nNetNode" };
    int dimid = 0;
    err = nc_inq_dimid(ncidp, mesh2dNodes.c_str(), &dimid);

    std::size_t num_nodes;
    auto read_name = new char[NC_MAX_NAME];
    err = nc_inq_dim(ncidp, dimid, read_name, &num_nodes);


    std::string mesh2dEdges{ "nNetLink" };
    int dimidedges = 0;
    err = nc_inq_dimid(ncidp, mesh2dEdges.c_str(), &dimid);

    std::size_t num_edges;
    err = nc_inq_dim(ncidp, dimid, read_name, &num_edges);

    std::vector<double> nodeX(num_nodes, 0.0);
    std::vector<double> nodeY(num_nodes, 0.0);
    std::vector<double> nodeZ(num_nodes, 0.0);
    std::string mesh2dNodeX{ "NetNode_x" };
    int varid = 0;
    err = nc_inq_varid(ncidp, mesh2dNodeX.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeX[0]);

    std::string mesh2dNodeY{ "NetNode_y" };
    err = nc_inq_varid(ncidp, mesh2dNodeY.c_str(), &varid);
    err = nc_get_var_double(ncidp, varid, &nodeY[0]);


    std::string mesh2dEdgeNodes{ "NetLink" };
    err = nc_inq_varid(ncidp, mesh2dEdgeNodes.c_str(), &varid);

    std::vector<int> edge_nodes(num_edges * 2, 0);
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
        edge_nodes[index] -= 1;
        edges[i].first = edge_nodes[index];
        
        index++;
        edge_nodes[index] -= 1;
        edges[i].second = edge_nodes[index] - 1;
        index++;
    }

    // now build node-edge mapping
    auto gridgeom = LoadLibrary("gridgeomStateful_dll.dll");
    auto start(std::chrono::steady_clock::now());

    int gridStateId;
    int ierr = ggeo_new_grid(gridStateId);
    if (ierr != 0) 
    {
        std::cout << "ggeo_new_grid failed " << std::endl;
        return 0;
    }

    MeshGeometryDimensions meshGeometryDimensions;
    meshGeometryDimensions.numnode = num_nodes;
    meshGeometryDimensions.numedge = num_edges;

    MeshGeometry meshGeometry;
    meshGeometry.nodex = &nodeX[0];
    meshGeometry.nodey = &nodeY[0];
    meshGeometry.nodez = &nodeZ[0];
    meshGeometry.edge_nodes = &edge_nodes[0];

    ierr = ggeo_set_state(gridStateId, meshGeometryDimensions, meshGeometry, false);

    if (ierr != 0)
    {
        std::cout << "ggeo_set_state failed " << std::endl;
        return 0;
    }

    int isTriangulationRequired = 0;
    int isAccountingForLandBoundariesRequired = 0;
    int projectToLandBoundaryOption = 0;
    OrthogonalizationParametersNative orthogonalizationParametersNative;
    GeometryListNative geometryListNativePolygon;
    GeometryListNative geometryListNativeLandBoundaries;
    ierr = ggeo_orthogonalize(gridStateId, isTriangulationRequired, isAccountingForLandBoundariesRequired, projectToLandBoundaryOption,
        orthogonalizationParametersNative, geometryListNativePolygon, geometryListNativeLandBoundaries);

    if (ierr != 0)
    {
        std::cout << "ggeo_orthogonalize failed " << std::endl;
        return 0;
    }

    auto end(std::chrono::steady_clock::now());
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();

    std::cout << "Elapsed time " << duration << " s " << std::endl;
    std::cout << "Test finished " << std::endl;

	return 0;
}
