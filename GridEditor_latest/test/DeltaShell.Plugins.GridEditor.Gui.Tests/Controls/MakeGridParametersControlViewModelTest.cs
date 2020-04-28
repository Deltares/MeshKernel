using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Gui.Controls;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controls
{
    [TestFixture]
    [Category(TestCategory.DataAccess)]
    public class MakeGridParametersControlViewModelTest
    {

        [Test]
        public void GivenMakeGridParametersControlViewModel_ReadAsciiFileDimensionsCommandWithoutValidPath_ShouldNotReadDimensions()
        {
            //Arrange
            var viewModel = new MakeGridParametersControlViewModel();

            // Act
            viewModel.AsciiFilePath = "abc";

            // Assert
            Assert.AreEqual("", viewModel.FileParametersString);
        }
    }
}