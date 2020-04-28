using System;
using System.Diagnostics;
using System.Threading;
using System.Windows;
using DelftTools.Controls;
using DelftTools.Shell.Gui;
using DelftTools.TestUtils;
using DelftTools.Utils.Reflection;
using DeltaShell.Plugins.GridEditor.Gui.Controllers;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.SharpMapGis.Gui.Forms;
using NUnit.Framework;
using Rhino.Mocks;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests
{
    [TestFixture]
    public class GridEditorGuiPluginTest
    {
        [Test]
        public void TestDefaultSettings()
        {
            var plugin = new GridEditorGuiPlugin();

            Assert.IsTrue(plugin?.Name?.Length > 0, "Plugin name should be specified");
            Assert.IsTrue(plugin?.DisplayName?.Length > 0, "Plugin display name should be specified");
            Assert.IsTrue(plugin?.Description?.Length > 0, "Plugin description should be specified");

            var versionInfo = FileVersionInfo.GetVersionInfo(typeof(GridEditorGuiPlugin).Assembly.Location);
            Assert.AreEqual(versionInfo.ProductVersion, plugin.Version);
            Assert.AreEqual(versionInfo.FileVersion, plugin.FileFormatVersion);
        }

        [Test]
        public void GraphicsProviderIsCreatedOnce()
        {
            var graphics = "/DeltaShell.Plugins.GridEditor.Gui;component/Controls/Dictionaries/IconBrushes.xaml";
            Application.LoadComponent(new Uri(graphics, UriKind.Relative));

            var plugin = new GridEditorGuiPlugin();

            var graphicsProvider = plugin.GraphicsProvider;

            Assert.AreEqual(graphicsProvider, plugin.GraphicsProvider, "GraphicsProvider should be instantiated once");
        }

        [Test]
        public void WhenChangingActiveView_CallsRibbonStopGridEditing()
        {
            var mocks = new MockRepository();
            var gui = mocks.StrictMock<IGui>();
            var viewList = mocks.StrictMock<IViewList>();
            var viewListToolWindows = mocks.StrictMock<IViewList>();

            var mapView = mocks.StrictMock<MapView>();
            var ribbon = mocks.StrictMock<IGridEditorRibbon>();

            gui.Expect(g => g.DocumentViews).Return(viewList).Repeat.Any();

            viewList.Expect(l => l.ActiveView).Return(mapView).Repeat.Any();
            viewListToolWindows.Expect(l => l.ActiveView).Return(null);

            ribbon.Expect(r => r.StopGridEditing());
            ribbon.Expect(r => r.ResetUnstructuredGrids());

            gui.Expect(g => g.ToolWindowViews).Return(viewListToolWindows).Repeat.Any();

            mocks.ReplayAll();

            var plugin = new GridEditorGuiPlugin { Gui = gui };

            TypeUtils.SetField(plugin, "gridEditorRibbon", ribbon);

            plugin.OnActiveViewChanged(mapView);

            plugin.Activate();

            plugin.OnActiveViewChanged(mapView);

            mocks.VerifyAll();
        }

        [Test]
        [Apartment(ApartmentState.STA)]
        public void RemovingMapViewCallsRibbonStopGridEditing()
        {
            var mocks = new MockRepository();
            var gui = mocks.StrictMock<IGui>();
            var viewList = mocks.StrictMock<IViewList>();
            var ribbon = mocks.StrictMock<IGridEditorRibbon>();

            var mapView = mocks.StrictMock<MapView>();
            var otherView = mocks.StrictMock<IView>();

            gui.Expect(g => g.DocumentViews).Return(viewList).Repeat.Any();
            
            viewList.Expect(l => l.ActiveView).Return(mapView).Repeat.Any();
            
            ribbon.Expect(r => r.StopGridEditing());

            mocks.ReplayAll();

            var plugin = new GridEditorGuiPlugin { Gui = gui };

            TypeUtils.SetField(plugin, "gridEditorRibbon", ribbon);
            
            // Should not trigger StopGridEditing
            plugin.OnViewRemoved(mapView);

            plugin.Activate();

            // Should not trigger StopGridEditing
            plugin.OnViewRemoved(otherView);

            // Should trigger StopGridEditing
            plugin.OnViewRemoved(mapView);

            mocks.VerifyAll();
        }

        [Test, Category(TestCategory.Integration)]
        [Apartment(ApartmentState.STA)]
        public void GivenGridEditorGuiPlugin_CallingRibbonCommandHandler_ShouldCreateRibbonWithMapInteractor()
        {
            //Arrange
            var mapview = new MapView();

            var mocks = new MockRepository();
            var gui = mocks.StrictMock<IGui>();
            var viewList = mocks.StrictMock<IViewList>();

            gui.Expect(g => g.DocumentViews).Return(viewList).Repeat.Any();
            viewList.Expect(l => l.ActiveView).Return(mapview).Repeat.Any();
            
            mocks.ReplayAll();

            var plugin = new GridEditorGuiPlugin { Gui = gui };
            
            // Act
            var ribbon = plugin.RibbonCommandHandler;
            Assert.IsInstanceOf<GridEditorRibbon>(ribbon);

            var mapcontrol = ((MapInteractor)((GridEditorRibbon)ribbon).MapInteractor).GetCurrentMapControl();

            // Assert
            Assert.AreEqual(mapview.MapControl, mapcontrol);
        }
    }
}