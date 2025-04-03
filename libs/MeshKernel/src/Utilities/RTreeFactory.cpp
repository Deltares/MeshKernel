#include "MeshKernel/Utilities/RTreeFactory.hpp"
#include "MeshKernel/Utilities/RTreeSphericalToCartesian.hpp"

std::unique_ptr<meshkernel::RTreeBase> meshkernel::RTreeFactory::Create(Projection projection)
{
    switch (projection)
    {
        using enum Projection;
    case cartesian:
        return std::make_unique<RTree<bg::cs::cartesian>>();
    case spherical:
    case sphericalAccurate:
        return std::make_unique<RTreeSphericalToCartesian>();
        // return std::make_unique<RTree<bg::cs::geographic<bg::degree>>>();
    default:
        throw MeshKernelError("Invalid projection '{}'", ProjectionToString(projection));
    }
}
