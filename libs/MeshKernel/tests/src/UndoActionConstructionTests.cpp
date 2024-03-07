#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/CompoundUndoAction.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/DeleteNodeAction.hpp"
#include "MeshKernel/ResetEdgeAction.hpp"
#include "MeshKernel/ResetNodeAction.hpp"
#include "MeshKernel/SphericalCoordinatesOffsetAction.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

#include "MockUndoAction.hpp"

namespace mk = meshkernel;

// Tests in this file test only basic functionality of the specific undo action
// Such as construction and that members are set correctly

TEST(UndoActionConstructionTests, AddEdgeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt edgeId = 10;
    mk::UInt start = 20;
    mk::UInt end = 30;

    std::unique_ptr<mk::AddEdgeAction> action = mk::AddEdgeAction::Create(*mesh, edgeId, start, end);
    EXPECT_EQ(edgeId, action->EdgeId());
    EXPECT_EQ(start, action->GetEdge().first);
    EXPECT_EQ(end, action->GetEdge().second);
}

TEST(UndoActionConstructionTests, DeleteEdgeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt edgeId = 10;
    mk::UInt start = 20;
    mk::UInt end = 30;

    std::unique_ptr<mk::DeleteEdgeAction> action = mk::DeleteEdgeAction::Create(*mesh, edgeId, start, end);
    EXPECT_EQ(edgeId, action->EdgeId());
    EXPECT_EQ(start, action->GetEdge().first);
    EXPECT_EQ(end, action->GetEdge().second);
}

TEST(UndoActionConstructionTests, ResetEdgeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt edgeId = 10;
    mk::UInt initStart = 20;
    mk::UInt initEnd = 30;
    mk::UInt updateStart = 21;
    mk::UInt updateEnd = 31;

    std::unique_ptr<mk::ResetEdgeAction> action = mk::ResetEdgeAction::Create(*mesh, edgeId, {initStart, initEnd}, {updateStart, updateEnd});
    EXPECT_EQ(edgeId, action->EdgeId());
    EXPECT_EQ(initStart, action->InitialEdge().first);
    EXPECT_EQ(initEnd, action->InitialEdge().second);
    EXPECT_EQ(updateStart, action->UpdatedEdge().first);
    EXPECT_EQ(updateEnd, action->UpdatedEdge().second);
}

TEST(UndoActionConstructionTests, AddNodeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt nodeId = 10;
    mk::Point node(0.5, 0.5);

    std::unique_ptr<mk::AddNodeAction> action = mk::AddNodeAction::Create(*mesh, nodeId, node);
    EXPECT_EQ(nodeId, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
}

TEST(UndoActionConstructionTests, DeleteNodeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt nodeId = 10;
    mk::Point node(0.5, 0.5);

    std::unique_ptr<mk::DeleteNodeAction> action = mk::DeleteNodeAction::Create(*mesh, nodeId, node);
    EXPECT_EQ(nodeId, action->NodeId());
    EXPECT_EQ(node.x, action->Node().x);
    EXPECT_EQ(node.y, action->Node().y);
    // TODO test the delete edge actions
}

TEST(UndoActionConstructionTests, ResetNodeTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    mk::UInt nodeId = 10;
    mk::Point initial(1.0, 2.0);
    mk::Point updated(2.0, 1.0);

    std::unique_ptr<mk::ResetNodeAction> action = mk::ResetNodeAction::Create(*mesh, nodeId, initial, updated);
    EXPECT_EQ(nodeId, action->NodeId());
    EXPECT_EQ(initial.x, action->InitialNode().x);
    EXPECT_EQ(initial.y, action->InitialNode().y);
    EXPECT_EQ(updated.x, action->UpdatedNode().x);
    EXPECT_EQ(updated.y, action->UpdatedNode().y);
}

TEST(UndoActionConstructionTests, SphericalCoordinatesOffsetActionTest)
{
    auto mesh = MakeRectangularMeshForTesting(2, 2, 1.0, meshkernel::Projection::cartesian);

    const double minX = 1.0;
    const double maxX = 5.0;

    std::unique_ptr<mk::SphericalCoordinatesOffsetAction> action = mk::SphericalCoordinatesOffsetAction::Create(*mesh, minX, maxX);

    std::vector<mk::UInt> nodesToIncrease{0, 1, 2, 4};
    std::vector<mk::UInt> nodesToDecrease{5, 6, 8};

    std::vector<mk::Point> originalNodes{{0.0, 0.0}, {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0}, {13.0, 14.0}, {15.0, 16.0}, {17.0, 18.0}};
    std::vector<mk::Point> transformedNodes(originalNodes);
    // The expected nodeal values after the offset has been applied
    std::vector<mk::Point> expectedNodes(originalNodes);

    // Generated the expected data.
    for (mk::UInt i = 0; i < nodesToIncrease.size(); ++i)
    {
        expectedNodes[nodesToIncrease[i]].x += 360.0;
    }

    for (mk::UInt i = 0; i < nodesToDecrease.size(); ++i)
    {
        expectedNodes[nodesToDecrease[i]].x -= 360.0;
    }

    // Set up offset
    for (mk::UInt nodeId : nodesToIncrease)
    {
        action->AddIncrease(nodeId);
    }

    for (mk::UInt nodeId : nodesToDecrease)
    {
        action->AddDecrease(nodeId);
    }

    EXPECT_EQ(action->MinX(), minX);
    EXPECT_EQ(action->MaxX(), maxX);

    constexpr double tolerance = 1.0e-10;

    action->ApplyOffset(transformedNodes);

    // Number of nodes should not have changed
    ASSERT_EQ(transformedNodes.size(), expectedNodes.size());

    for (mk::UInt i = 0; i < expectedNodes.size(); ++i)
    {
        EXPECT_NEAR(transformedNodes[i].x, expectedNodes[i].x, tolerance);
        EXPECT_NEAR(transformedNodes[i].y, expectedNodes[i].y, tolerance);
    }

    action->UndoOffset(transformedNodes);

    // Number of nodes should not have changed
    ASSERT_EQ(transformedNodes.size(), originalNodes.size());

    for (mk::UInt i = 0; i < originalNodes.size(); ++i)
    {
        EXPECT_NEAR(transformedNodes[i].x, originalNodes[i].x, tolerance);
        EXPECT_NEAR(transformedNodes[i].y, originalNodes[i].y, tolerance);
    }
}

class TestUndoAction : public mk::UndoAction
{
public:
    /// @brief How many objects were created
    static int s_created;
    /// @brief How many objects are committed
    static int s_committed;
    /// @brief How many objects are restored
    static int s_restored;

    static std::unique_ptr<TestUndoAction> Create(const int value, std::vector<int>& undoActionValues)
    {
        return std::make_unique<TestUndoAction>(value, undoActionValues);
    }

    TestUndoAction(const int value, std::vector<int>& undoActionValues) : m_value(value), m_undoActionValues(undoActionValues)
    {
        m_undoActionValues[s_created] = m_value;
        ++s_created;
        // Objects must be created in the committed state.
        ++s_committed;
    }

    int Value() const
    {
        return m_value;
    }

private:
    void DoCommit() override
    {
        m_undoActionValues[s_committed] = m_value;
        ++s_committed;
        --s_restored;
    }

    void DoRestore() override
    {
        m_undoActionValues[s_restored] = m_value;
        --s_committed;
        ++s_restored;
    }

    int m_value;
    std::vector<int>& m_undoActionValues;
};

int TestUndoAction::s_created = 0;
int TestUndoAction::s_committed = 0;
int TestUndoAction::s_restored = 0;

TEST(UndoActionConstructionTests, CompoundUndoActionTest)
{
    std::unique_ptr<mk::CompoundUndoAction> compoundAction = mk::CompoundUndoAction::Create();

    const std::vector<int> expectedValues{1, 10, 3, 21};
    // Values will be assigned in the order depending on the creation, restore or committed call
    // creation and committed should be in same order as expectedValues.
    // restored values should be in reverse.
    std::vector<int> undoActionValues(expectedValues.size());

    for (int value : expectedValues)
    {
        compoundAction->Add(TestUndoAction::Create(value, undoActionValues));
    }

    // Expect values to be in created order
    for (size_t i = 0; i < undoActionValues.size(); ++i)
    {
        EXPECT_EQ(undoActionValues[i], expectedValues[i]);
    }

    ////////////////////////////////
    // Check action has been constructed correctly
    ////////////////////////////////
    ASSERT_EQ(std::distance(compoundAction->begin(), compoundAction->end()),
              static_cast<std::iterator_traits<mk::CompoundUndoAction::const_iterator>::difference_type>(expectedValues.size()));

    ASSERT_EQ(TestUndoAction::s_created, static_cast<int>(expectedValues.size()));
    EXPECT_EQ(TestUndoAction::s_committed, static_cast<int>(expectedValues.size()));
    EXPECT_EQ(TestUndoAction::s_restored, 0);

    size_t count = 0;

    for (auto iter = compoundAction->begin(); iter != compoundAction->end(); ++iter)
    {
        const TestUndoAction* testAction = dynamic_cast<const TestUndoAction*>(iter->get());
        EXPECT_EQ(testAction->Value(), expectedValues[count]);
        ++count;
    }

    ////////////////////////////////
    // Check restore action is performed correctly
    ////////////////////////////////
    compoundAction->Restore();

    // Expect values to be in reverse created order
    for (size_t i = 0; i < undoActionValues.size(); ++i)
    {
        EXPECT_EQ(undoActionValues[i], expectedValues[undoActionValues.size() - i - 1]);
    }

    // Number of created TestUndoAction object should not change
    EXPECT_EQ(TestUndoAction::s_created, static_cast<int>(expectedValues.size()));
    // After restore, there should be zero committed TestUndoAction objects
    EXPECT_EQ(TestUndoAction::s_committed, 0);
    // After restore all TestUndoAction objects should be restored
    EXPECT_EQ(TestUndoAction::s_restored, static_cast<int>(expectedValues.size()));

    ////////////////////////////////
    // Check commit action is performed correctly
    ////////////////////////////////
    compoundAction->Commit();

    // Expect values to be in created order
    for (size_t i = 0; i < undoActionValues.size(); ++i)
    {
        EXPECT_EQ(undoActionValues[i], expectedValues[i]);
    }

    // Number of created TestUndoAction object should not change
    ASSERT_EQ(TestUndoAction::s_created, static_cast<int>(expectedValues.size()));
    // After commit all TestUndoAction objects should be committed
    EXPECT_EQ(TestUndoAction::s_committed, static_cast<int>(expectedValues.size()));
    // After restore, there should be zero restored TestUndoAction objects
    EXPECT_EQ(TestUndoAction::s_restored, 0);
}

TEST(UndoActionConstructionTests, CompoundUndoActionAddRestoredTest)
{
    std::unique_ptr<mk::CompoundUndoAction> compoundAction = mk::CompoundUndoAction::Create();

    std::unique_ptr<mk::UndoAction> simpleAction = std::make_unique<MockUndoAction>();
    simpleAction->Restore();

    EXPECT_THROW(compoundAction->Add(std::move(simpleAction)), mk::ConstraintError);
}
