#include "MeshKernel/RemoveDisconnectedRegions.hpp"
#include "MeshKernel/Exceptions.hpp"

#include <queue>

meshkernel::UInt meshkernel::RemoveDisconnectedRegions::GetNeighbour(const std::array<UInt, 2>& edge, const UInt elementId) const
{
    if (elementId == constants::missing::uintValue)
    {
        return constants::missing::uintValue;
    }

    if (edge[0] == elementId)
    {
        return edge[1];
    }

    if (edge[1] == elementId)
    {
        return edge[0];
    }

    return constants::missing::uintValue;
}

void meshkernel::RemoveDisconnectedRegions::LabelConnectedRegion(const Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, const UInt unlabledElementId, UInt& elementCount) const
{
    // Use a flood fill like algorithm to label all elements, connected via edges, with the same region id.

    std::queue<UInt> toBeProcessed;
    // Add a single unlabeled element to the queue
    toBeProcessed.push(unlabledElementId);
    elementRegionId[unlabledElementId] = regionId;
    elementCount = 1;

    // Process the queue, until there are no nodes left to process.
    while (!toBeProcessed.empty())
    {
        const UInt currentElement = toBeProcessed.front();
        toBeProcessed.pop();

        // Get the neighbours of this node from the element-edge connectivity graph.
        // Only across faces
        for (const UInt edge : mesh.m_facesEdges[currentElement])
        {
            const UInt neighbour = GetNeighbour(mesh.m_edgesFaces[edge], currentElement);

            if (neighbour != constants::missing::uintValue && elementRegionId[neighbour] == constants::missing::uintValue)
            {
                // If any neighbour element is labeled undefined it will get current region id
                // and added to the end of the queue to be processed later.
                elementRegionId[neighbour] = regionId;
                ++elementCount;
                toBeProcessed.push(neighbour);
            }
        }
    }
}

void meshkernel::RemoveDisconnectedRegions::LabelSingleDomainRegion(const Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, UInt& elementCount) const
{
    UInt unlabledElementId = constants::missing::uintValue;
    elementCount = 0;

    // Find an element that has not yet been labeled, any as yet unlabeled element will do,
    // so use the first one that we find in the list.
    for (UInt i = 0; i < elementRegionId.size(); ++i)
    {
        // Find first unlabeled element
        if (elementRegionId[i] == constants::missing::uintValue)
        {
            unlabledElementId = i;
            break;
        }
    }

    // If such an element is found then label all connected elements (across faces)
    // with the current regionId.
    if (unlabledElementId != constants::missing::uintValue)
    {
        LabelConnectedRegion(mesh, regionId, elementRegionId, unlabledElementId, elementCount);
    }
}

void meshkernel::RemoveDisconnectedRegions::LabelAllDomainRegions(const Mesh2D& mesh, std::vector<UInt>& elementRegionId, std::vector<std::pair<UInt, UInt>>& regionCount) const
{
    bool allRegionsLabeled = false;
    // Initialise with the first region identifier.
    UInt regionId = 1;

    elementRegionId.resize(mesh.GetNumFaces(), constants::missing::uintValue);
    regionCount.clear();

    while (!allRegionsLabeled)
    {
        UInt elementCount = 0;
        LabelSingleDomainRegion(mesh, regionId, elementRegionId, elementCount);

        if (elementCount > 0)
        {
            // If more than 0 elements were labaled, then save this region-id and number of elements.
            regionCount.emplace_back(regionId, elementCount);
            ++regionId;
        }
        else
        {
            // If no elements were found then we are complete.
            allRegionsLabeled = true;
        }
    }
}

void meshkernel::RemoveDisconnectedRegions::RemoveDetachedRegions(Mesh2D& mesh, const UInt regionId, std::vector<UInt>& elementRegionId, UInt& numberOfElementsRemoved) const
{
    numberOfElementsRemoved = 0;

    // Loop over all element region ids, deleting all those that do not have the same region id as the master id.
    for (UInt i = 0; i < elementRegionId.size(); ++i)
    {
        if (elementRegionId[i] != regionId)
        {
            ++numberOfElementsRemoved;
            std::ranges::for_each(mesh.m_facesEdges[i], [&mesh](const UInt edge)
                                  { mesh.DeleteEdge(edge); });
            elementRegionId[i] = constants::missing::uintValue;
        }
    }
}

void meshkernel::RemoveDisconnectedRegions::Compute(Mesh2D& mesh) const
{
    // Label the elements of each discontiguous region of the mesh with a unique identifier for the region.
    // A mapping between the identifier for the region and the number of elements in the same region is retained
    // Only the region containing the largest number of elements will be kept, all other regions removed.

    // Map from element id (the domain, index of the array) to region id.
    std::vector<UInt> elementRegionId;
    // List of number of elements labeled with number
    // the domain (index) is the region id and the range, regionCount[id], is the element count.
    std::vector<std::pair<UInt, UInt>> regionCount;

    // Generate region identifier for all elements in the mesh and a mapping between the region id the the number of elements.
    LabelAllDomainRegions(mesh, elementRegionId, regionCount);

    // If more that 1 region has been found then remove all regions that have fewer elements
    // than the region with the maximum number of elements.
    if (regionCount.size() > 1)
    {
        // We assume the region with the largest number of elements is the region of interest we would like to retain.
        // all other regions can (and will) be removed.
        UInt maxRegionElementCount = 0;
        UInt mainRegionId = 0;

        for (const auto& [regionId, regionElementCount] : regionCount)
        {
            if (regionElementCount > maxRegionElementCount)
            {
                maxRegionElementCount = regionElementCount;
                mainRegionId = regionId;
            }
        }

        UInt numberOfElementsRemoved = 0;
        // Remove all elements from the regions that do not have the main region id.
        RemoveDetachedRegions(mesh, mainRegionId, elementRegionId, numberOfElementsRemoved);
        mesh.Administrate();
    }
}
