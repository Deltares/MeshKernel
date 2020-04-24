#pragma once

namespace GridGeomApi
{
    struct MeshGeometryDimensions
    {
        char name[255];
        int dim;
        int numnode;
        int numedge;
        int numface;
        int maxnumfacenodes;
        int numlayer;
        int layertype;
        int nnodes;
        int nbranches;
        int ngeometry;
        int epgs;
    };
}


