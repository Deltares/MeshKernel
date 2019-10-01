#pragma once

namespace GridGeomApi
{
struct GeometryListNative
{
    int type;
    double geometrySeperator;
    double innerOuterSeperator;
    int numberOfCoordinates;
    double* xCoordinates = nullptr;
    double* yCoordinates = nullptr;
    double* zCoordinates = nullptr;
};
}


