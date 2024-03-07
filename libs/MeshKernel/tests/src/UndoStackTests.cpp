#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <iterator>

#include "MeshKernel/AddEdgeAction.hpp"
#include "MeshKernel/AddNodeAction.hpp"
#include "MeshKernel/DeleteEdgeAction.hpp"
#include "MeshKernel/DeleteNodeAction.hpp"
#include "MeshKernel/ResetEdgeAction.hpp"
#include "MeshKernel/ResetNodeAction.hpp"
#include "MeshKernel/UndoAction.hpp"
#include "MeshKernel/UndoActionStack.hpp"

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
    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(0);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);

    // Should call the DoRestore for the last action added to the stack, which is undoAction2
    undoActionStack.Undo();
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoCommit()).Times(1);

    // Should call the DoCommit for the last action that was undone
    undoActionStack.Commit();
    EXPECT_EQ(undoActionStack.Size(), 2);

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);
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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoCommit()).Times(0);

    // There should be no actions in the restored state.
    EXPECT_FALSE(undoActionStack.Commit());
    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction4->State(), mk::UndoAction::Committed);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction4->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Undo());
    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction4->State(), mk::UndoAction::Restored);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

    EXPECT_CALL(*rawUndoAction1, DoRestore()).Times(0);
    EXPECT_CALL(*rawUndoAction2, DoRestore()).Times(1);
    EXPECT_CALL(*rawUndoAction3, DoRestore()).Times(1);

    // Should call the DoRestore for the last two actions added to the stack
    undoActionStack.Undo();
    undoActionStack.Undo();
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);
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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);

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

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Restored);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Commit()); // action2.Restore
    EXPECT_EQ(undoActionStack.Size(), 3);

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Restored);

    EXPECT_TRUE(undoActionStack.Undo());  // action2.Undo
    EXPECT_TRUE(undoActionStack.Undo());  // action1.Undo
    EXPECT_FALSE(undoActionStack.Undo()); // no action undone
    EXPECT_EQ(undoActionStack.Size(), 3);

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

    EXPECT_TRUE(undoActionStack.Undo());    // action3.Undo
    EXPECT_TRUE(undoActionStack.Undo());    // action2.Undo
    EXPECT_TRUE(undoActionStack.Commit());  // action2.Restore
    EXPECT_TRUE(undoActionStack.Commit());  // action3.Restore
    EXPECT_FALSE(undoActionStack.Commit()); // no action restored

    EXPECT_EQ(rawUndoAction1->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction2->State(), mk::UndoAction::Committed);
    EXPECT_EQ(rawUndoAction3->State(), mk::UndoAction::Committed);
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
    for (mk::UInt i = 0; i < mk::UndoActionStack::MaxUndoSize; ++i)
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
    for (mk::UInt i = 0; i < mk::UndoActionStack::MaxUndoSize; ++i)
    {
        undoActionStack.Add(std::make_unique<MockUndoAction>());
        EXPECT_EQ(undoActionStack.Size(), i + 1);
    }

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), mk::UndoActionStack::MaxUndoSize);

    // Now add another undo action. The number on the stack should not increase.
    undoActionStack.Add(std::make_unique<MockUndoAction>());

    // Check the size
    EXPECT_EQ(undoActionStack.Size(), mk::UndoActionStack::MaxUndoSize);

    // Clear all undo actions form stack
    undoActionStack.Clear();
    EXPECT_EQ(undoActionStack.Size(), 0);
}
