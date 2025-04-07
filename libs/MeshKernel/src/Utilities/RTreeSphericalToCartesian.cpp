#include "MeshKernel/Utilities/RTreeSphericalToCartesian.hpp"

#include <cmath>

meshkernel::RTreeSphericalToCartesian::Point2D meshkernel::RTreeSphericalToCartesian::convert(const Point& node) const
{
    double x = constants::geometric::earth_radius * std::cos(node.y * constants::conversion::degToRad) * std::cos(node.x * constants::conversion::degToRad);
    double y = constants::geometric::earth_radius * std::cos(node.y * constants::conversion::degToRad) * std::sin(node.x * constants::conversion::degToRad);
    double z = constants::geometric::earth_radius * std::sin(node.y * constants::conversion::degToRad);
    return {x, y, z};
}

void meshkernel::RTreeSphericalToCartesian::Search(Point const& node, double searchRadiusSquared, bool findNearest)
{
    if (Empty())
    {
        throw AlgorithmError("RTree is empty, search cannot be performed");
    }

    m_queryCache.reserve(m_queryVectorCapacity);
    m_queryCache.clear();

    const Point2D nodeSought = convert(node);
    const auto searchRadius = std::sqrt(searchRadiusSquared);
    Box2D const box(Point2D(bg::get<0>(nodeSought) - searchRadius, bg::get<1>(nodeSought) - searchRadius, bg::get<2>(nodeSought) - searchRadius),
                    Point2D(bg::get<0>(nodeSought) + searchRadius, bg::get<1>(nodeSought) + searchRadius, bg::get<2>(nodeSought) + searchRadius));

    auto pointIsNearby = [&nodeSought, &searchRadiusSquared](Value2D const& v)
    { return bg::comparable_distance(v.first, nodeSought) <= searchRadiusSquared; };

    if (findNearest)
    {
        m_rtree2D.query(bgi::within(box) && bgi::satisfies(pointIsNearby) && bgi::nearest(nodeSought, 1),
                        std::back_inserter(m_queryCache));
    }
    else
    {
        m_rtree2D.query(bgi::within(box) && bgi::satisfies(pointIsNearby),
                        std::back_inserter(m_queryCache));
    }

    m_queryIndices.clear();

    if (findNearest && !m_queryCache.empty())
    {
        m_queryIndices.emplace_back(m_queryCache[0].second);
    }
    else
    {
        for (const auto& entry : m_queryCache)
        {
            m_queryIndices.emplace_back(entry.second);
        }
    }
}

void meshkernel::RTreeSphericalToCartesian::SearchPoints(Point const& node, double searchRadiusSquared)
{
    Search(node, searchRadiusSquared, false);
}

void meshkernel::RTreeSphericalToCartesian::SearchNearestPoint(Point const& node, double searchRadiusSquared)
{
    Search(node, searchRadiusSquared, true);
}

void meshkernel::RTreeSphericalToCartesian::SearchNearestPoint(Point const& node)
{
    if (Empty())
    {
        throw AlgorithmError("RTree is empty, search cannot be performed");
    }

    m_queryCache.reserve(m_queryVectorCapacity);
    m_queryCache.clear();

    const Point2D nodeSought = convert(node);
    m_rtree2D.query(bgi::nearest(nodeSought, 1), std::back_inserter(m_queryCache));

    if (!m_queryCache.empty())
    {
        m_queryIndices.clear();
        m_queryIndices.emplace_back(m_queryCache[0].second);
    }
}

void meshkernel::RTreeSphericalToCartesian::DeleteNode(UInt position)
{
    if (Empty())
    {
        throw AlgorithmError("RTree is empty, deletion cannot performed");
    }

    if (const auto numberRemoved = m_rtree2D.remove(m_points[position]); numberRemoved != 1)
    {
        return;
    }
    m_points[position] = {Point2D{constants::missing::doubleValue, constants::missing::doubleValue, constants::missing::doubleValue}, std::numeric_limits<UInt>::max()};
}
