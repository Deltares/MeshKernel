using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controls
{
    [TestFixture]
    [Category(TestCategory.DataAccess)]
    public class MakeGridParametersControlViewModelTest
    {
        [Test]
        public void GivenMakeGridParametersControlViewModel_ReadAsciiFileDimensionsCommand_ShouldReadDimensions()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("SchematisatieInt9x9.asc");
            var viewModel = new MakeGridParametersControlViewModel
            {
                GetFilePath = ()=> path
            };

            // Act
            viewModel.ReadAsciiFileDimensionsCommand.Execute(null);

            // Assert
            var expectedFileParametersString = "Type = Square\r\n" + 
                                               "BlockSize = 20\r\n" + 
                                               "X BlockSize = 20\r\n" +
                                               "Y BlockSize = 20\r\n" + 
                                               "Origin X Coordinate = 0\r\n" +
                                               "Origin Y Coordinate = 0\r\n" + 
                                               "Number Of Rows  = 8\r\n" +
                                               "Number Of Columns = 8\r\n";

            Assert.AreEqual(expectedFileParametersString, viewModel.FileParametersString);
        }

        [Test]
        public void GivenMakeGridParametersControlViewModel_ReadAsciiFileDimensionsCommandWithForcePointsOnCellCenters_ShouldReadDimensions()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("SchematisatieInt9x9.asc");
            var viewModel = new MakeGridParametersControlViewModel
            {
                GetFilePath = () => path,
                ForcePointsOnCellCenters = true
            };

            // Act
            viewModel.ReadAsciiFileDimensionsCommand.Execute(null);

            // Assert
            var expectedFileParametersString = "Type = Square\r\n" +
                                               "BlockSize = 20\r\n" +
                                               "X BlockSize = 20\r\n" +
                                               "Y BlockSize = 20\r\n" +
                                               "Origin X Coordinate = -10\r\n" +
                                               "Origin Y Coordinate = -10\r\n" +
                                               "Number Of Rows  = 9\r\n" +
                                               "Number Of Columns = 9\r\n";

            Assert.AreEqual(expectedFileParametersString, viewModel.FileParametersString);
        }

        [Test]
        public void GivenMakeGridParametersControlViewModel_ReadAsciiFileDimensionsCommandWithoutValidPath_ShouldNotReadDimensions()
        {
            //Arrange
            var viewModel = new MakeGridParametersControlViewModel
            {
                GetFilePath = () => "abc"
            };

            // Act
            viewModel.ReadAsciiFileDimensionsCommand.Execute(null);

            // Assert
            Assert.AreEqual("", viewModel.FileParametersString);
        }


        [Test]
        public void GivenMakeGridParametersControlViewModel_SetDimensionsFromAsciiCommand_ShouldSetReadDimensionsToGridParameters()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("SchematisatieInt9x9.asc");
            var viewModel = new MakeGridParametersControlViewModel
            {
                GridParameters = MakeGridParameters.CreateDefault(),
                AsciiFilePath = path,
                ForcePointsOnCellCenters = true
            };

            var parameters = viewModel.GridParameters;

            Assert.AreEqual(parameters.XGridBlockSize, 10);
            Assert.AreEqual(parameters.YGridBlockSize, 10);
            Assert.AreEqual(parameters.NumberOfColumns, 3);
            Assert.AreEqual(parameters.NumberOfRows, 3);

            // Act
            Assert.IsTrue(viewModel.SetDimensionsFromAsciiCommand.CanExecute(null));

            viewModel.SetDimensionsFromAsciiCommand.Execute(null);

            // Assert
            Assert.AreEqual(parameters.XGridBlockSize, 0);
            Assert.AreEqual(parameters.YGridBlockSize, 0);
            Assert.AreEqual(parameters.NumberOfColumns, 0);
            Assert.AreEqual(parameters.NumberOfRows, 0);
        }
    }
}