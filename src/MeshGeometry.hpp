#pragma once

namespace GridGeomApi
{
    struct MeshGeometry
    {
        int* edge_nodes = nullptr;
        int* face_nodes = nullptr;
        int* edge_faces = nullptr;
        int* face_edges = nullptr;
        int* face_links = nullptr;

        double* nnodex = nullptr;
        double* nnodey = nullptr;
        int* nedge_nodes = nullptr;
        double* nbranchlengths = nullptr;
        int* nbranchgeometrynodes = nullptr;
        double* ngeopointx = nullptr;
        double* ngeopointy = nullptr;
        int* nbranchorder = nullptr;
        int* branchidx = nullptr;
        double* branchoffsets = nullptr;

        double* nodex = nullptr;
        double* nodey = nullptr;
        double* nodez = nullptr;
        double* edgex = nullptr;
        double* edgey = nullptr;
        double* edgez = nullptr;
        double* facex = nullptr;
        double* facey = nullptr;
        double* facez = nullptr;

        double* layer_zs = nullptr;
        double* interface_zs = nullptr;
        int     startIndex = 0;
    };
}


