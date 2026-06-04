#include "MeshKernel/TriangulationGenerator.hpp"

#include <algorithm>
#include <execution>
#include <tuple>

std::unique_ptr<meshkernel::Mesh2D> meshkernel::SimpleTriangulationGenerator::generate(const Polygons& polygon) const
{

    if (polygon.GetNumPolygons() != 1)
    {
        throw MeshKernelError("Cannot generate a triangulation for {} polygons", polygon.GetNumPolygons());
    }

    // generate samples in the first polygonal enclosure
    auto const generatedPoints = polygon.Enclosure(0).GeneratePoints(scaleFactor_ < 0.0 ? meshkernel::constants::missing::doubleValue : scaleFactor_);

    return std::make_unique<Mesh2D>(generatedPoints, polygon, polygon.GetProjection());
}

double meshkernel::SepranTriangulationGenerator::minimumEdgeDelta(const std::vector<meshkernel::Point>& polygonNodes)
{
    double delta = std::numeric_limits<double>::max();

    for (size_t i = 0; i + 1 < polygonNodes.size(); ++i)
    {
        double dx = polygonNodes[i + 1].x - polygonNodes[i].x;
        double dy = polygonNodes[i + 1].y - polygonNodes[i].y;

        delta = std::min(delta, std::sqrt(dx * dx + dy * dy));
    }

    return delta;
}

std::vector<meshkernel::Point> meshkernel::SepranTriangulationGenerator::pointsFromFlatArray(const std::vector<double>& coordinates, const int numberOfPoints)
{
    std::vector<meshkernel::Point> meshNodes(numberOfPoints);

    for (int i = 0; i < numberOfPoints; ++i)
    {
        meshNodes[i].x = coordinates[2 * i];
        meshNodes[i].y = coordinates[2 * i + 1];
    }

    return meshNodes;
}

std::vector<std::reference_wrapper<const meshkernel::Polygon>> meshkernel::SepranTriangulationGenerator::generatePolygonReferences(const Polygons& polygon)
{
    const auto& enclosure = polygon.Enclosure(0);

    std::vector<std::reference_wrapper<const meshkernel::Polygon>> boundaryLoops;
    boundaryLoops.reserve(1 + enclosure.NumberOfInner());
    boundaryLoops.emplace_back(enclosure.Outer());

    for (meshkernel::UInt i = 0; i < enclosure.NumberOfInner(); ++i)
    {
        boundaryLoops.emplace_back(enclosure.Inner(i));
    }

    return boundaryLoops;
}

std::tuple<std::vector<meshkernel::Edge>, std::vector<std::vector<meshkernel::UInt>>, std::vector<std::uint8_t>>
meshkernel::SepranTriangulationGenerator::gatherEdgesAndFaces(const std::vector<int>& kmeshc, const int numberOfElements)
{

    std::vector<meshkernel::Edge> edges(3 * numberOfElements);
    std::vector<std::vector<meshkernel::UInt>> faceNodes(numberOfElements);
    std::vector<std::uint8_t> numFaceNodes(numberOfElements, 0);

    // Gather all edges in the mesh
    // The first stage will find shared edges twice, once for each triangle. The duplicate will be removed later.
    for (int i = 0; i < numberOfElements; ++i)
    {
        int index = 3 * i;

        meshkernel::UInt n1 = static_cast<meshkernel::UInt>(kmeshc[index] - 1);
        meshkernel::UInt n2 = static_cast<meshkernel::UInt>(kmeshc[index + 1] - 1);
        meshkernel::UInt n3 = static_cast<meshkernel::UInt>(kmeshc[index + 2] - 1);

        // Which is better, reserve and push back or allocate and assign?
        edges[index] = n1 < n2 ? meshkernel::Edge(n1, n2) : meshkernel::Edge(n2, n1);
        edges[index + 1] = n2 < n3 ? meshkernel::Edge(n2, n3) : meshkernel::Edge(n3, n2);
        edges[index + 2] = n3 < n1 ? meshkernel::Edge(n3, n1) : meshkernel::Edge(n1, n3);

        // Set the nodes of the element.
        faceNodes[i].resize(3);
        faceNodes[i][0] = n1;
        faceNodes[i][1] = n2;
        faceNodes[i][2] = n3;

        // Indicate that there are 3 nodes for this element
        numFaceNodes[i] = 3;
    }

    // Remove the duplicated edges.
    // std::pair has a predefined less-than (Edge is a std;:pair<UInt>)
    std::sort(std::execution::par, edges.begin(), edges.end());
    auto [first, last] = std::ranges::unique(edges);
    edges.erase(first, last);

    return {edges, faceNodes, numFaceNodes};
}

std::unique_ptr<meshkernel::Mesh2D> meshkernel::SepranTriangulationGenerator::generate(const Polygons& polygon) const
{

    if (polygon.GetNumPolygons() != 1)
    {
        throw MeshKernelError("Cannot generate a triangulation for {} polygons", polygon.GetNumPolygons());
    }

    std::vector<std::reference_wrapper<const Polygon>> boundaryLoops(generatePolygonReferences(polygon));

    const int numberOfBoundaryNodes = std::accumulate(boundaryLoops.begin(), boundaryLoops.end(), 0, [](int sum, const auto& poly)
                                                      { return sum + static_cast<int>(poly.get().Size()); });

    const int dimension = 2;
    const int elementIdentifier = 3; // 3-node linear triangles
    const int numberOfPolygons = static_cast<int>(boundaryLoops.size());
    const int numberOfHoles = numberOfPolygons - 1;

    // 1. Interleave coordinates into bcord: [x1, y1, x2, y2...]
    // Flatten the outer polygon followed by each inner polygon, without separators.
    std::vector<double> boundaryCoordinates(2 * numberOfBoundaryNodes);

    // 2. Map edgeNodeConnectivity: Fortran uses 1-based indexing
    std::vector<int> edgeNodeConnectivity(numberOfBoundaryNodes, 0);

    // 3. Map boundary segments: Fortran boundary(2, numberOfPolygons)
    // One curve entry per closed loop, flattened in column-major order.

    std::vector<int> boundaryConnectivity(2 * numberOfPolygons);

    // 4. Compute elementSizing array for local polyline point matching
    const int numberElementSizing = 0;         // numberOfBoundaryNodes;
    std::vector<double> elementSizing(1, 0.0); // 3 * numberOfBoundaryNodes);

    int pointOffset = 0;

    for (int loopIndex = 0; loopIndex < numberOfPolygons; ++loopIndex)
    {
        const auto& loop = boundaryLoops[loopIndex].get().Nodes();

        boundaryConnectivity[2 * loopIndex] = loopIndex + 1;
        boundaryConnectivity[2 * loopIndex + 1] = pointOffset + 1;

        for (int i = 0; i < static_cast<int>(loop.size()); ++i)
        {
            boundaryCoordinates[2 * (pointOffset + i)] = loop[i].x;
            boundaryCoordinates[2 * (pointOffset + i) + 1] = loop[i].y;
            edgeNodeConnectivity[pointOffset + i] = pointOffset + i + 1;
        }

        pointOffset += static_cast<int>(loop.size());
    }

    // 5. Dynamic Memory Allocation for Output Buffers
    // Get estimate of area covered by polygon, subtracting area covered by holes
    const double estimatedArea = polygon.Enclosure(0).ComputeSurfaceArea();
    // Compute the minimum spacing between points of the polygon, both outer and all inner polygons
    const auto [minimumDelta, maximumDeltaDummy] = polygon.Enclosure(0).SegmentLengthExtrema();

    const int estimatedNumberOfElements = static_cast<int>((estimatedArea / (0.433 * minimumDelta * minimumDelta)) * 3.5);
    const int maximumNumberOfNodes = numberOfBoundaryNodes + 3 * estimatedNumberOfElements;
    const int maximumNumberOfElements = 2 * maximumNumberOfNodes;

    std::vector<double> triangulationNodes(dimension * maximumNumberOfNodes, 0.0);
    std::vector<int> triangulationElementNodes(elementIdentifier * maximumNumberOfElements, 0);

    // 6. Dummies and Placeholders
    std::vector<int> holeinfo(2 * (numberOfHoles + 2), 0);
    int forcedControlPoints = 0;
    int surfaceSequenceNumber = 1;
    int auxiliaryAlignment = 0;

    std::vector<int> numnodextcurvs(1, 0);
    std::vector<int> curvenumbers(1, 0);
    std::vector<double> rinput(1, 0.0);
    std::vector<int> userpoints(1, 0);

    // 7. Make the Call
    int newMesh = 1;
    int numberOfPoints = maximumNumberOfNodes;
    int numberOfElements = maximumNumberOfElements;

    mshoce_(&newMesh, triangulationNodes.data(), triangulationElementNodes.data(), &elementIdentifier, &numberOfBoundaryNodes, boundaryCoordinates.data(),
            edgeNodeConnectivity.data(), boundaryConnectivity.data(), &numberOfPolygons, &numberOfPoints, &numberOfElements,
            holeinfo.data(), &numberOfHoles, &numberElementSizing, elementSizing.data(), userpoints.data(),
            &surfaceSequenceNumber, &auxiliaryAlignment, numnodextcurvs.data(), curvenumbers.data(),
            rinput.data(), &forcedControlPoints, &dimension);

    // Recover array of Points
    std::vector<Point> meshNodes(pointsFromFlatArray(triangulationNodes, numberOfPoints));

    // Recover arrays of edges, face-node connectivity and number of nodes per face.
    auto [edges, faceNodes, numFaceNodes] = gatherEdgesAndFaces(triangulationElementNodes, numberOfElements);

    return std::make_unique<Mesh2D>(edges, meshNodes, faceNodes, numFaceNodes, polygon.GetProjection());
}
