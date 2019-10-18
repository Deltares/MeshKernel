using System;
using System.Collections.Generic;
using System.Linq;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Helpers;
using NetTopologySuite.Extensions.Grids;
using NUnit.Framework;
using Rhino.Mocks;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controls.Ribbon
{
    [TestFixture]
    public class GridEditorRibbonViewModelTest
    {
        [Test]
        public void CheckShowNodes()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.ShowVertices = true),
                r => r.ShowVertices = true, 
                nameof(GridEditorRibbonViewModel.ShowVertices));
        }

        [Test]
        public void CheckNodeSize()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.VertexSize = 5),
                r => r.VertexSize = 5, 
                nameof(GridEditorRibbonViewModel.VertexSize));
        }

        [Test]
        public void CheckDrawFilledInPolygons()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.DrawFilledInPolygons = true),
                r => r.DrawFilledInPolygons = true,
                nameof(GridEditorRibbonViewModel.DrawFilledInPolygons));
        }

        [Test]
        public void CheckPolygonRefinementDistance()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.PolygonRefinementDistance = 15),
                r => r.PolygonRefinementDistance = 15,
                nameof(GridEditorRibbonViewModel.PolygonRefinementDistance));
        }

        [Test]
        public void CheckPolygonOffsetDistance()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.PolygonOffsetDistance = 5),
                r => r.PolygonOffsetDistance = 5,
                nameof(GridEditorRibbonViewModel.PolygonOffsetDistance));
        }

        [Test]
        public void CheckDrawSelectedVertices()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.DrawSelectedVertices = true),
                r => r.DrawSelectedVertices = true,
                nameof(GridEditorRibbonViewModel.DrawSelectedVertices));
        }

        [Test]
        public void CallingRefresh_ShouldGeneratePropertyChanged()
        {
            //Arrange
            var expectedPropertyNames = new[]
            {
                nameof(GridEditorRibbonViewModel.PolygonToolEnabled),
                nameof(GridEditorRibbonViewModel.ChangePolygonToolEnabled),
                nameof(GridEditorRibbonViewModel.SplineToolEnabled),
                nameof(GridEditorRibbonViewModel.LandBoundariesToolEnabled),
                nameof(GridEditorRibbonViewModel.InsertEdgesToolEnabled),

                nameof(GridEditorRibbonViewModel.ProjectToLandBoundaryOption),
                nameof(GridEditorRibbonViewModel.IsAccountingForLandBoundariesRequired),
                nameof(GridEditorRibbonViewModel.IsTriangulationRequired),

                nameof(GridEditorRibbonViewModel.SplinesToCurvilinearParameters),
                nameof(GridEditorRibbonViewModel.CurvilinearParameters),
                nameof(GridEditorRibbonViewModel.PolygonRefinementDistance),
                nameof(GridEditorRibbonViewModel.MakeGridParameters),
                nameof(GridEditorRibbonViewModel.OrthogonalizationParameters),
                nameof(GridEditorRibbonViewModel.DeleteMeshOption),

                nameof(GridEditorRibbonViewModel.VertexSize),
                nameof(GridEditorRibbonViewModel.DrawSelectedVertices),
                nameof(GridEditorRibbonViewModel.DrawFilledInPolygons),
                nameof(GridEditorRibbonViewModel.ShowVertices),

                nameof(GridEditorRibbonViewModel.IsEditing),
            };

            var viewModel = new GridEditorRibbonViewModel();

            // Act
            var propertyChangedCount = 0;
            viewModel.PropertyChanged += (s, a) =>
            {
                Assert.Contains(a.PropertyName, expectedPropertyNames);
                propertyChangedCount++;
            };

            viewModel.RefreshState();

            // Assert
            Assert.AreEqual(propertyChangedCount, expectedPropertyNames.Length);
        }

        [Test]
        public void CheckIsEditing()
        {
            CheckPropertyChangedAndControllerSet(
                controller =>
                {
                    controller.Expect(c => c.IsEditing = true);
                    controller.Expect(c => c.IsEditing).Return(false);
                },
                r =>
                {
                    r.IsEditing = true;
                }, 
                nameof(GridEditorRibbonViewModel.IsEditing));
        }

        [Test]
        public void SettingIsTriangulationRequired_ShouldTriggerControllerIsTriangulationRequired()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.IsTriangulationRequired = true),
                r => r.IsTriangulationRequired = true,
                nameof(GridEditorRibbonViewModel.IsTriangulationRequired));
        }
        
        [Test]
        public void SettingIsAccountingForLandBoundariesRequired_ShouldTriggerControllerIsAccountingForLandBoundariesRequired()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.IsAccountingForLandBoundariesRequired = true),
                r => r.IsAccountingForLandBoundariesRequired = true,
                nameof(GridEditorRibbonViewModel.IsAccountingForLandBoundariesRequired));
        }

        [Test]
        public void SettingProjectToLandBoundaryOption_ShouldTriggerControllerProjectToLandBoundaryOption()
        {
            CheckPropertyChangedAndControllerSet(
                controller => controller.Expect(c => c.ProjectToLandBoundaryOption = ProjectToLandBoundaryOptions.No),
                r => r.ProjectToLandBoundaryOption = ProjectToLandBoundaryOptions.No,
                nameof(GridEditorRibbonViewModel.ProjectToLandBoundaryOption));
        }

        [Test]
        public void IsEditingShouldNotCreatePropertyChangedWhenTheValueIsTheSame()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.IsEditing).Return(true).Repeat.Any();

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            var propertyChangedCount = 0;
            ribbonViewModel.PropertyChanged += (sender, args) => { propertyChangedCount++; };

            ribbonViewModel.IsEditing = true;

            Assert.AreEqual(0, propertyChangedCount, "Setting the same value should not create a property changed event");

            mocks.VerifyAll();
        }

        private static void CheckPropertyChangedAndControllerSet(Action<IGridEditorController> setControllerExpectation,
            Action<GridEditorRibbonViewModel> viewModelAction, string propertyName)
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            setControllerExpectation(controller);
            
            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            var propertyChangedCount = 0;
            ribbonViewModel.PropertyChanged += (sender, args) =>
            {
                propertyChangedCount++;
                Assert.AreEqual(args.PropertyName, propertyName);
            };

            viewModelAction(ribbonViewModel);

            Assert.AreEqual(1, propertyChangedCount, $"Changing property should create property changed event for {propertyName}");

            mocks.VerifyAll();
        }

        [Test]
        public void SettingSelectedGridCreatesTwoPropertyChangedEvents()
        {
            var ribbonViewModel = new GridEditorRibbonViewModel();

            var propertyChangedCount = 0;
            ribbonViewModel.PropertyChanged += (sender, args) =>
            {
                propertyChangedCount++;
                Assert.IsTrue(new[]
                {
                    nameof(GridEditorRibbonViewModel.SelectedGrid),
                    nameof(GridEditorRibbonViewModel.CanEdit)
                }.Contains(args.PropertyName), "Should be expected property");
            };

            Assert.IsFalse(ribbonViewModel.CanEdit);

            ribbonViewModel.SelectedGrid = new UnstructuredGrid();

            Assert.IsTrue(ribbonViewModel.CanEdit);

            Assert.AreEqual(2, propertyChangedCount, $"Changing property should create two property changed events");
        }

        [Test]
        public void SettingUnstructuredGridsSetsFirstGridAsSelectedGrid()
        {
            var ribbonViewModel = new GridEditorRibbonViewModel();

            var propertyChangedCount = 0;
            ribbonViewModel.PropertyChanged += (sender, args) =>
            {
                propertyChangedCount++;
                Assert.IsTrue(new[]
                {
                    nameof(GridEditorRibbonViewModel.UnstructuredGrids),
                    nameof(GridEditorRibbonViewModel.SelectedGrid),
                    nameof(GridEditorRibbonViewModel.CanEdit),
                    nameof(GridEditorRibbonViewModel.MoreThanOneGrid)
                }.Contains(args.PropertyName), "Should be expected property");
            };

            Assert.IsNull(ribbonViewModel.SelectedGrid, "Selected grid is initially not set");

            ribbonViewModel.UnstructuredGrids = new List<UnstructuredGrid>
            {
                new UnstructuredGrid(),
                new UnstructuredGrid()
            };

            Assert.AreEqual(ribbonViewModel.SelectedGrid, ribbonViewModel.UnstructuredGrids[0], "Selected grid is initially not set");

            Assert.AreEqual(4, propertyChangedCount, $"Changing property should create three property changed events");
        }

        [Test]
        public void DeleteVerticesCommandShouldTriggerControllerDeleteSelectedVertices()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.DeleteSelectedVertices());

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            ribbonViewModel.DeleteVerticesCommand.Execute(null);

            mocks.VerifyAll();
        }

        [Test]
        public void FlipEdgesCommandShouldTriggerControllerFlipEdges()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.FlipEdges());

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            ribbonViewModel.FlipEdgesCommand.Execute(null);

            mocks.VerifyAll();
        }

        [Test]
        public void MergeTwoVerticesCommandShouldTriggerControllerMergeVertices()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.MergeTwoVertices());

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            ribbonViewModel.MergeTwoVerticesCommand.Execute(null);

            mocks.VerifyAll();
        }

        [Test]
        public void MergeVerticesCommandShouldTriggerControllerMergeVertices()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.MergeVertices());

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            ribbonViewModel.MergeVerticesCommand.Execute(null);

            mocks.VerifyAll();
        }

        [Test]
        [TestCase(ImportExportType.Grid)]
        [TestCase(ImportExportType.Polygons)]
        [TestCase(ImportExportType.Spline)]
        [TestCase(ImportExportType.LandBoundary)]
        [TestCase(ImportExportType.Samples)]
        public void GridImportCommandShouldTriggerControllerImportExportWithCorrectType(ImportExportType importExportType)
        {
            var path = "test.txt";

            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.ImportExport(ImportExportAction.Import, importExportType, path));

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller)
            {
                GetFilePath = (s, action) => path
            };

            ribbonViewModel.ImportCommand.Execute(importExportType);

            mocks.VerifyAll();
        }

        [Test]
        [TestCase(ImportExportType.Grid)]
        [TestCase(ImportExportType.Polygons)]
        [TestCase(ImportExportType.Spline)]
        [TestCase(ImportExportType.LandBoundary)]
        [TestCase(ImportExportType.Samples)]
        public void GridExportCommandShouldTriggerControllerImportExportWithCorrectType(ImportExportType importExportType)
        {
            var path = "test.txt";

            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.ImportExport(ImportExportAction.Export, importExportType, path));

            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller)
            {
                GetFilePath = (s, action) => path
            };

            ribbonViewModel.ExportCommand.Execute(importExportType);

            mocks.VerifyAll();
        }

        [Test]
        public void GridExportCommandShouldNotTriggerControllerImportExportIfPathIsEmpty()
        {
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();
            
            mocks.ReplayAll();

            var ribbonViewModel = new GridEditorRibbonViewModel(controller);

            ribbonViewModel.ExportCommand.Execute(ImportExportType.Grid);

            mocks.VerifyAll();
        }

        [Test]
        public void CallingSwitchWorkFlowCommand_ShouldCauseMultiplePropertyChangedEventsAndSetSelectedGenerateGridWorkFlowType()
        {
            //Arrange
            var propertyChangedCount = 0;
            var viewModel = new GridEditorRibbonViewModel();
            viewModel.PropertyChanged += (s, a) =>
            {
                propertyChangedCount++;
                Assert.Contains(a.PropertyName, new []
                {
                    nameof(viewModel.SelectedGenerateGridWorkFlowType),
                    nameof(viewModel.SelectedGenerateGridWorkFlowTypeString),
                    nameof(viewModel.IsPolygonToolVisible),
                    nameof(viewModel.IsSplineToolVisible)
                });
            };

            // Act
            viewModel.SwitchWorkFlowCommand.Execute(GenerateGridWorkFlowType.Regular);

            // Assert
            Assert.AreEqual(4, propertyChangedCount);
            Assert.AreEqual(GenerateGridWorkFlowType.Regular, viewModel.SelectedGenerateGridWorkFlowType);
        }
    }
}