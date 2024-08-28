#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iterator>
#include <random>

#include "MeshKernel/UndoActions/AddEdgeAction.hpp"
#include "MeshKernel/UndoActions/AddNodeAction.hpp"
#include "MeshKernel/UndoActions/DeleteEdgeAction.hpp"
#include "MeshKernel/UndoActions/DeleteNodeAction.hpp"
#include "MeshKernel/UndoActions/ResetEdgeAction.hpp"
#include "MeshKernel/UndoActions/ResetNodeAction.hpp"
#include "MeshKernel/UndoActions/UndoAction.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"

#include "MeshKernel/Constants.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Mesh2D.hpp"

#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"

#include "MockUndoAction.hpp"

namespace mk = meshkernel;

TEST(UndoStackTests, CheckSimpleUndo)
{
    // A test on the basic undo functionality

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    EXPECT_EQ(undoActionStack.Size(), 2);

    // Check state is correct
    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(0);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);

    // Still 2 items
    EXPECT_EQ(undoActionStack.Size(), 2);
}

TEST(UndoStackTests, CheckSimpleUndoWithRedo)
{
    // A test on the basic undo and commit (redo) functionality

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);

    // Should call the DoCommit for the last action that was undone
    undoActionStack.Commit();
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
}

TEST(UndoStackTests, CheckMultipleUndo)
{
    // A test on the basic multiple undo's functionality

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(undoActionStack.Size(), 3);
}

TEST(UndoStackTests, CheckMultipleAllUndo)
{
    // A test on the basic multiple undo's functionality
    // And no further undo is performed if requested

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(0);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_TRUE(undoActionStack.Undo());
    // Should not perform any undo action
    EXPECT_FALSE(undoActionStack.Undo());

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);
}

TEST(UndoStackTests, CheckMultipleUndoWithSingleRedo)
{
    // A test on the basic multiple undo's with single redo

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(0);

    // Should call the DoCommit for the last action that was undone
    undoActionStack.Commit();

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);
}

TEST(UndoStackTests, CheckMultipleUndoWithSingleRedoIntermediateAdd)
{
    // A test on the basic multiple undo's with single redo

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();
    auto undoAction4 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction4 = undoAction4.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction4, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();
    undoActionStack.Undo();
    // At this point, the raw pointers rawUndoAction2 and rawUndoAction3 are invlaid.
    // Since the unique_ptr's have been destroyed and, therefore, contents deallocated.
    rawUndoAction2 = nullptr;
    rawUndoAction3 = nullptr;

    undoActionStack.Add(std::move(undoAction4));
    // The two undo actions that have been restored should have been removed
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);

    // There should be no actions in the restored state.
    EXPECT_FALSE(undoActionStack.Commit());
    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction4->GetState(), mk::UndoAction::State::Committed);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction4->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction4->GetState(), mk::UndoAction::State::Restored);

    // Attempt to undo.
    // There should be no other actions to undo.
    EXPECT_FALSE(undoActionStack.Undo());
    EXPECT_EQ(undoActionStack.Size(), 2);
}

TEST(UndoStackTests, CheckMultipleUndoWithMultipleRedo)
{
    // A test on the basic multiple undo's with multiple redo's
    // Should return to the original state with all actions committed.

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    undoActionStack.Undo();
    undoActionStack.Undo();
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(1);

    // Should call the DoCommit for the last two action that was undone
    EXPECT_TRUE(undoActionStack.Commit());
    EXPECT_EQ(undoActionStack.Size(), 3);
    EXPECT_TRUE(undoActionStack.Commit());
    EXPECT_EQ(undoActionStack.Size(), 3);
    // Should not perform any other commit (redo) actions
    EXPECT_FALSE(undoActionStack.Commit());
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);
}

TEST(UndoStackTests, CheckMultipleUndoCycles)
{
    // A test on the basic multiple undo's with multiple redo's cycles
    // Should return to the original state with all actions committed.

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Have to use raw pointers, since the pointer value will be "moved" from the unique_ptr when adding to the undo stack
    MockUndoAction* rawUndoAction1 = undoAction1.get();
    auto undoAction2 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction2 = undoAction2.get();
    auto undoAction3 = std::make_unique<MockUndoAction>();
    MockUndoAction* rawUndoAction3 = undoAction3.get();

    undoActionStack.Add(std::move(undoAction1));
    undoActionStack.Add(std::move(undoAction2));
    undoActionStack.Add(std::move(undoAction3));
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(3);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(2);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(1);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(3);
    EXPECT_CALL(*rawUndoAction3, DoCommit()).Times(2);

    // Should call the DoRestore for the last two actions added to the stack
    EXPECT_TRUE(undoActionStack.Undo()); // action3.Undo
    EXPECT_TRUE(undoActionStack.Undo()); // action2.Undo
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action2.Restore
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Undo());  // action2.Undo
    EXPECT_TRUE(undoActionStack.Undo());  // action1.Undo
    EXPECT_FALSE(undoActionStack.Undo()); // no action undone
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action1.Restore

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Restored);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action2.Restore

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action3.Restore

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);

    EXPECT_TRUE(undoActionStack.Undo());    // action3.Undo
    EXPECT_TRUE(undoActionStack.Undo());    // action2.Undo
    EXPECT_TRUE(undoActionStack.Commit());  // action2.Restore
    EXPECT_TRUE(undoActionStack.Commit());  // action3.Restore
    EXPECT_FALSE(undoActionStack.Commit()); // no action restored

    EXPECT_EQ(rawUndoAction1->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction2->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(rawUndoAction3->GetState(), mk::UndoAction::State::Committed);
    EXPECT_EQ(undoActionStack.Size(), 3);
}

TEST(UndoStackTests, AddingRestoredAction)
{

    mk::UndoActionStack undoActionStack;

    auto undoAction1 = std::make_unique<MockUndoAction>();
    // Restore action before adding to stack
    undoAction1->Restore();

    EXPECT_THROW(undoActionStack.Add(std::move(undoAction1)), mk::ConstraintError);
}

TEST(UndoStackTests, MemorySizeStackTest)
{
    mk::UndoActionStack undoActionStack;
    std::uint64_t expectedSize = sizeof(undoActionStack);

    // Fill undo actions to maximum
    for (mk::UInt i = 0; i < mk::UndoActionStack::DefaultMaxUndoSize; ++i)
    {
        auto undoAction = std::make_unique<MockUndoAction>();
        expectedSize += undoAction->MemorySize();
        undoActionStack.Add(std::move(undoAction));
    }

    // Check the memory-size
    EXPECT_EQ(undoActionStack.MemorySize(), expectedSize);

    // Clear all undo actions form stack
    undoActionStack.Clear();
    EXPECT_EQ(undoActionStack.Size(), 0);
}

TEST(UndoStackTests, ExeedingMaximumUndoActions)
{
    mk::UndoActionStack undoActionStack;

    EXPECT_EQ(undoActionStack.Size(), 0);

    // Fill undo actions to maximum
    for (mk::UInt i = 0; i < mk::UndoActionStack::DefaultMaxUndoSize; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>());
        EXPECT_EQ(undoActionStack.Size(), i + 1);
    }

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), mk::UndoActionStack::DefaultMaxUndoSize);

    // Now add another undo action. The number on the stack should not increase.
    undoActionStack.Add(std::make_unique<MockUndoAction>());

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), mk::UndoActionStack::DefaultMaxUndoSize);

    // Clear all undo actions form stack
    undoActionStack.Clear();
    EXPECT_EQ(undoActionStack.Size(), 0);
}

TEST(UndoStackTests, RemovingUndoActions)
{
    mk::UndoActionStack undoActionStack;
    int actionId1 = 1;
    int actionId2 = 2;
    int actionId3 = 3;
    int unknownActionId = 4;

    mk::UInt actionCount = (mk::UndoActionStack::DefaultMaxUndoSize - mk::UndoActionStack::DefaultMaxUndoSize % 3) / 3;

    EXPECT_EQ(undoActionStack.Size(), 0);

    // Fill undo actions
    for (mk::UInt i = 0; i < actionCount; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>(), actionId1);
        undoActionStack.Add(std::make_unique<MockUndoAction>(), actionId2);
        undoActionStack.Add(std::make_unique<MockUndoAction>(), actionId3);
    }

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), 3 * actionCount);

    // Remove all undo actions associated with actionId2
    EXPECT_EQ(undoActionStack.Remove(actionId2), actionCount);

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), 2 * actionCount);

    // There are no undo actions associated with unknownActionId. So, should remove 0 undo actions
    EXPECT_EQ(undoActionStack.Remove(unknownActionId), 0);

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), 2 * actionCount);

    // Clear all undo actions form stack
    undoActionStack.Clear();
    EXPECT_EQ(undoActionStack.Size(), 0);
}

TEST(UndoStackTests, RemovingUndoActionsRandomised)
{
    constexpr mk::UInt NumberOfUndoActions = 4;

    mk::UndoActionStack undoActionStack;
    std::array<int, NumberOfUndoActions> actionIds{1, 10, 21, 33};
    std::array<int, NumberOfUndoActions> actionCounts{0, 0, 0, 0};
    int unknownActionId = 45;

    mk::UInt actionCount = (mk::UndoActionStack::DefaultMaxUndoSize - mk::UndoActionStack::DefaultMaxUndoSize % NumberOfUndoActions) / NumberOfUndoActions;
    mk::UInt totalActionCount = NumberOfUndoActions * actionCount;

    EXPECT_EQ(undoActionStack.Size(), 0);

    std::mt19937 generator;

    // Fill undo actions
    for (mk::UInt i = 0; i < totalActionCount; ++i)
    {
        mk::UInt randomNum = static_cast<mk::UInt>(generator()) % NumberOfUndoActions;
        undoActionStack.Add(std::make_unique<MockUndoAction>(), actionIds[randomNum]);
        ++actionCounts[randomNum];
    }

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), NumberOfUndoActions * actionCount);

    // Remove all undo actions associated with actionId[1]
    EXPECT_EQ(undoActionStack.Remove(actionIds[1]), actionCounts[1]);
    totalActionCount -= actionCounts[1];

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), totalActionCount);

    // Remove all undo actions associated with actionId[3]
    EXPECT_EQ(undoActionStack.Remove(actionIds[3]), actionCounts[3]);
    totalActionCount -= actionCounts[3];

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), totalActionCount);

    // There are no undo actions associated with unknownActionId. So, should remove 0 undo actions
    EXPECT_EQ(undoActionStack.Remove(unknownActionId), 0);

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), totalActionCount);

    // Clear all undo actions form stack
    undoActionStack.Clear();
    EXPECT_EQ(undoActionStack.Size(), 0);
}

TEST(UndoStackTests, InitaliseStackSizeToZero)
{
    mk::UndoActionStack undoActionStack(0);
    EXPECT_EQ(undoActionStack.Size(), 0);

    // Add large number of undo actions
    for (mk::UInt i = 0; i < mk::UndoActionStack::DefaultMaxUndoSize; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>());
    }

    EXPECT_EQ(undoActionStack.Size(), 0);
    // There should be nothing to undo
    EXPECT_FALSE(undoActionStack.Undo());
}

TEST(UndoStackTests, SetStackSizeToZero)
{
    mk::UndoActionStack undoActionStack;

    EXPECT_EQ(undoActionStack.Size(), 0);
    undoActionStack.SetMaximumSize(0);

    // Fill undo actions to maximum
    for (mk::UInt i = 0; i < mk::UndoActionStack::DefaultMaxUndoSize; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>());
    }

    EXPECT_EQ(undoActionStack.Size(), 0);

    // Cannot undo, there are no items stored in stack
    EXPECT_FALSE(undoActionStack.Undo());
}

TEST(UndoStackTests, SetStackSizeToZeroAfterAdding)
{
    mk::UInt InitialStackSize = 20;
    // This number must be small than InitialStackSize but larger than 0
    mk::UInt RevisedStackSize = 10;

    mk::UndoActionStack undoActionStack(InitialStackSize);

    // Fill undo actions to maximum
    for (mk::UInt i = 0; i < InitialStackSize; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>());
    }

    EXPECT_EQ(undoActionStack.Size(), InitialStackSize);
    undoActionStack.SetMaximumSize(RevisedStackSize);
    EXPECT_EQ(undoActionStack.Size(), RevisedStackSize);

    for (mk::UInt i = 0; i < RevisedStackSize; ++i)
    {
        EXPECT_TRUE(undoActionStack.Undo());
    }

    // Cannot undo further, there are no items left in stack
    EXPECT_FALSE(undoActionStack.Undo());

    // Redo all the undone actions
    for (mk::UInt i = 0; i < RevisedStackSize; ++i)
    {
        EXPECT_TRUE(undoActionStack.Commit());
    }

    EXPECT_EQ(undoActionStack.Size(), RevisedStackSize);

    // Now finally set the maximum size to be zero.
    undoActionStack.SetMaximumSize(0);
    EXPECT_EQ(undoActionStack.Size(), 0);

    // Cannot undo, there are no items in stack
    EXPECT_FALSE(undoActionStack.Undo());
}
