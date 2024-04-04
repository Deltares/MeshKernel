#include "MeshKernel/CurvilinearGrid/ResetCurvilinearNodeAction.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"

std::unique_ptr<meshkernel::ResetCurvilinearNodeAction> meshkernel::ResetCurvilinearNodeAction::Create(CurvilinearGrid& grid,
                                                                                                       const CurvilinearGridNodeIndices nodeId,
                                                                                                       const Point& initial,
                                                                                                       const Point& updated,
                                                                                                       const bool recalculateNodeTypes)
{
    return std::make_unique<ResetCurvilinearNodeAction>(grid, nodeId, initial, updated, recalculateNodeTypes);
}

meshkernel::ResetCurvilinearNodeAction::ResetCurvilinearNodeAction(CurvilinearGrid& grid,
                                                                   const CurvilinearGridNodeIndices nodeId,
                                                                   const Point& initial,
                                                                   const Point& updated,
                                                                   const bool recalculateNodeTypes) : BaseMeshUndoAction<ResetCurvilinearNodeAction, CurvilinearGrid>(grid),
                                                                                                      m_nodeId(nodeId),
                                                                                                      m_initialNode(initial),
                                                                                                      m_updatedNode(updated),
                                                                                                      m_recalculateNodeTypes(recalculateNodeTypes) {}

meshkernel::CurvilinearGridNodeIndices meshkernel::ResetCurvilinearNodeAction::NodeId() const
{
    return m_nodeId;
}

const meshkernel::Point& meshkernel::ResetCurvilinearNodeAction::InitialNode() const
{
    return m_initialNode;
}

const meshkernel::Point& meshkernel::ResetCurvilinearNodeAction::UpdatedNode() const
{
    return m_updatedNode;
}

bool meshkernel::ResetCurvilinearNodeAction::RecalculateNodeTypes() const
{
    return m_recalculateNodeTypes;
}

void meshkernel::ResetCurvilinearNodeAction::Print(std::ostream& out) const
{
    out << "ResetCurvilinearNodeAction: state " << to_string(GetState()) << ", nodeId = {" << m_nodeId.m_n << ", " << m_nodeId.m_m << "} "
        << ", initial = {" << m_initialNode.x << ", " << m_initialNode.y << "} "
        << ", updated = {" << m_updatedNode.x << ", " << m_updatedNode.y << "}"
        << ", node-types " << std::boolalpha << m_recalculateNodeTypes
        << std::endl;
}
