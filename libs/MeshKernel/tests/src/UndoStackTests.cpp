#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/DeleteNodeAction.hpp"
#include "MeshKernel/MoveNodeAction.hpp"
#include "MeshKernel/ResetEdgeAction.hpp"
#include "MeshKernel/ResetNodeAction.hpp"
#include "MeshKernel/UndoActionStack.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

namespace mk = meshkernel;
using ::testing::AtLeast;

class MockUndoAction : public mk::UndoAction
{
public:
    MOCK_METHOD(void, DoCommit, (), (override));

    MOCK_METHOD(void, DoRestore, (), (override));
};

TEST(UndoStackTests, CheckSimpleUndo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));

    // Check state is correct
    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
}

TEST(UndoStackTests, CheckSimpleUndoWithRedo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);

    // Should call the DoCommit for the last action that was undone
    undoActionStack.Commit();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
}

TEST(UndoStackTests, CheckMultipleUndo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    std::unique_ptr<MockUndoAction> undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);
}

TEST(UndoStackTests, CheckMultipleAllUndo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    std::unique_ptr<MockUndoAction> undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());
    // Should not perform any undo action
    EXPECT_FALSE(undoActionStack.Undo());

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);
}

TEST(UndoStackTests, CheckMultipleUndoWithSingleRedo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    std::unique_ptr<MockUndoAction> undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(0);

    // Should call the DoCommit for the last action that was undone
    undoActionStack.Commit();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);
}

TEST(UndoStackTests, CheckMultipleUndoWithMultipleRedo)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    std::unique_ptr<MockUndoAction> undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    undoActionStack.Undo();
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(1);

    // Should call the DoCommit for the last two action that was undone
    EXPECT_TRUE(undoActionStack.Commit());
    EXPECT_TRUE(undoActionStack.Commit());
    // Should not perform any other commit (redo) actions
    EXPECT_FALSE(undoActionStack.Commit());

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);
}

TEST(UndoStackTests, CheckMultipleUndoCycles)
{
    mk::UndoActionStack undoActionStack;

    std::unique_ptr<MockUndoAction> undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    std::unique_ptr<MockUndoAction> undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    std::unique_ptr<MockUndoAction> undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(2);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(2);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    EXPECT_TRUE(undoActionStack.Undo()); // action3.Undo
    EXPECT_TRUE(undoActionStack.Undo()); // action2.Undo

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action2.Restore

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Undo()); // action2.Undo
    EXPECT_TRUE(undoActionStack.Undo()); // action1.Undo

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action1.Restore

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action2.Restore

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action3.Restore

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);
}
