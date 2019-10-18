using System;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using NetTopologySuite.Extensions.Grids;
using NUnit.Framework;
using Rhino.Mocks;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controls.Ribbon
{
    [TestFixture]
    public class GridEditorRibbonTest
    {
        [Test]
        [STAThread]
        public void IsGridRibbonVisibleCheckShouldCheckTabNameAndActiveView()
        {
            var grid = new UnstructuredGrid();
            var mocks = new MockRepository();
            var interactor = mocks.StrictMock<IMapInteractor>();

            interactor.Expect(i => i.GetUnstructuredGrids()).Return(new[] {grid});

            mocks.ReplayAll();

            var ribbon = new GridEditorRibbon
            {
                MapInteractor = interactor
            };

            Assert.IsTrue(ribbon.IsContextualTabVisible("geospatialContextualGroup", "GridEditorTab"));
            Assert.IsFalse(ribbon.IsContextualTabVisible("OtherGroupName", "GridEditorTab"));
            Assert.IsFalse(ribbon.IsContextualTabVisible("geospatialContextualGroup", "OtherTabName"));

            mocks.VerifyAll();

            mocks.BackToRecord(interactor);
            interactor.Expect(i => i.GetUnstructuredGrids()).Return(new UnstructuredGrid[] { });
            interactor.Replay();

            Assert.IsFalse(ribbon.IsContextualTabVisible("geospatialContextualGroup", "GridEditorTab"));

            mocks.VerifyAll();
        }

        [Test]
        [STAThread]
        public void StopGridEditing_ShouldNotTrigger_GridEditorControllerIsEditingState_WhenNotInEditMode()
        {
            //Arrange
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.IsEditing).Return(false).Repeat.Any();
            
            mocks.ReplayAll();

            // Act
            var ribbonViewModel = new GridEditorRibbon
            {
                GridEditingGroupViewModel = new GridEditorRibbonViewModel(controller)
            };

            ribbonViewModel.StopGridEditing();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        [STAThread]
        public void StopGridEditing_ShouldTrigger_GridEditorControllerIsEditingState()
        {
            //Arrange
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.IsEditing).Return(true).Repeat.Any();
            controller.Expect(c => c.IsEditing = false);

            mocks.ReplayAll();

            // Act
            var ribbonViewModel = new GridEditorRibbon
            {
                GridEditingGroupViewModel = new GridEditorRibbonViewModel(controller)
            };

            ribbonViewModel.StopGridEditing();

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        [STAThread]
        public void ResetUnstructuredGrids_ShouldTrigger_GridEditorControllerResetUnstructuredGrids()
        {
            //Arrange
            var mocks = new MockRepository();
            var controller = mocks.StrictMock<IGridEditorController>();

            controller.Expect(c => c.IsEditing).Return(true).Repeat.Any();
            controller.Expect(c => c.ResetUnstructuredGrids());

            mocks.ReplayAll();

            // Act
            var ribbonViewModel = new GridEditorRibbon
            {
                GridEditingGroupViewModel = new GridEditorRibbonViewModel(controller)
            };

            ribbonViewModel.ResetUnstructuredGrids();

            // Assert
            mocks.VerifyAll();
        }
    }
}