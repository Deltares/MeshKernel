using DeltaShell.Plugins.GridEditor.Data;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.Data
{
    [TestFixture]
    public class GridEditorStateTest
    {
        [Test]
        public void GridEditorData_ShouldHaveDefaultValues_WhenCreated()
        {
            //Arrange
            var data = new GridEditorState();

            // Assert
            Assert.NotNull(data.LandBoundaries);
            Assert.NotNull(data.Splines);
            Assert.NotNull(data.OrthogonalizationParameters);
        }
    }
}