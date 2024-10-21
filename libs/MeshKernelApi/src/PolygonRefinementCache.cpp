#include "MeshKernelApi/PolygonRefinementCache.hpp"

meshkernelapi::PolygonRefinementCache::PolygonRefinementCache(const std::vector<meshkernel::Point>& polyPoints,
                                                              const int firstIndex,
                                                              const int secondIndex,
                                                              const double edgeLength,
                                                              const std::vector<meshkernel::Point>& refinedPoints)
    : CachedPointValues(refinedPoints),
      m_polygonPoints(polyPoints),
      m_firstNodeIndex(firstIndex),
      m_secondNodeIndex(secondIndex),
      m_targetEdgeLength(edgeLength) {}

bool meshkernelapi::PolygonRefinementCache::ValidOptions(const std::vector<meshkernel::Point>& polyPoints,
                                                         const int firstIndex,
                                                         const int secondIndex,
                                                         const double edgeLength) const
{
    return firstIndex == m_firstNodeIndex &&
           secondIndex == m_secondNodeIndex &&
           edgeLength == m_targetEdgeLength &&
           polyPoints.size() == m_polygonPoints.size() &&
           std::equal(polyPoints.begin(), polyPoints.end(), m_polygonPoints.begin());
}
