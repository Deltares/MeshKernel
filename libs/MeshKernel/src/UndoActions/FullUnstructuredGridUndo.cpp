#include "MeshKernel/UndoActions/FullUnstructuredGridUndo.hpp"
#include "MeshKernel/Mesh.hpp"

std::unique_ptr<meshkernel::FullUnstructuredGridUndo> meshkernel::FullUnstructuredGridUndo::Create(Mesh& mesh)
{
    return std::make_unique<FullUnstructuredGridUndo>(mesh);
}

meshkernel::FullUnstructuredGridUndo::FullUnstructuredGridUndo(Mesh& mesh) : BaseMeshUndoAction<FullUnstructuredGridUndo, Mesh>(mesh), m_savedNodes(mesh.Nodes()), m_savedEdges(mesh.Edges()) {}

void meshkernel::FullUnstructuredGridUndo::Swap(std::vector<Point>& nodes, std::vector<Edge>& edges)
{
    std::swap(nodes, m_savedNodes);
    std::swap(edges, m_savedEdges);
}

std::uint64_t meshkernel::FullUnstructuredGridUndo::MemorySize() const
{
    return sizeof(FullUnstructuredGridUndo) + sizeof(Point) * m_savedNodes.capacity() + sizeof(Edge) * m_savedEdges.capacity();
}
