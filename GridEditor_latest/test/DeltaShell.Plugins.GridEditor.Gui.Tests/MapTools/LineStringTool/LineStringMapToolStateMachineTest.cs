using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool
{
    [TestFixture]
    public class LineStringMapToolStateMachineTest
    {
        [Test]
        public void GivenPolygonMapToolStateMachine_switchingState_ShouldWork()
        {
            //Arrange
            var stateMachine = new LineStringMapToolStateMachine();

            // Act
            Assert.AreEqual(LineStringMapToolState.Selecting,stateMachine.CurrentState);

            stateMachine.MoveNext(LineStringMapToolStateCommand.MouseDownNotOnFeatureNoSelection);

            // Assert
            Assert.AreEqual(LineStringMapToolState.AddingNew, stateMachine.CurrentState);
        }
    }
}