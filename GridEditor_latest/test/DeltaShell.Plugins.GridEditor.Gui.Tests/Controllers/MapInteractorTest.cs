using System;
using System.Collections.Generic;
using System.Linq;
using DelftTools.Utils.Collections.Generic;
using DelftTools.Utils.Reflection;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.Controllers;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapTools;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Grids;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.Layers;
using SharpMap.UI.Forms;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.Controllers
{
    [TestFixture]
    public class MapInteractorTest
    {
        [Test]
        public void GivenMapInteractor_CallingStartStopInteraction_ShouldAddRemoveMapToolsAndLayers()
        {
            //Arrange
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var tools = mocks.StrictMock<IList<IMapTool>>();
            var map = mocks.StrictMock<IMap>();
            var layers = mocks.StrictMock<IEventedList<ILayer>>();

            mapControl.Expect(c => c.Map).Return(map).Repeat.Times(3);
            mapControl.Expect(c => c.Tools).Return(tools).Repeat.Times(4);

            tools.Expect(t => t.Add(null)).IgnoreArguments().WhenCalled(m =>
            {
                var isCorrectType = m.Arguments[0] is MoveTool ||
                                    m.Arguments[0] is IGridEditorMapTool ||
                                    m.Arguments[0] is InsertVerticesMapTool ||
                                    m.Arguments[0] is DeleteEdgesMapTool ||
                                    m.Arguments[0] is MoveVerticesMapTool;
                Assert.True(isCorrectType);
            }).Repeat.Times(7);

            tools.Expect(t => t.Remove(null)).IgnoreArguments().WhenCalled(m =>
            {
                var isCorrectType = m.Arguments[0] is MoveTool ||
                                    m.Arguments[0] is IGridEditorMapTool;
                Assert.True(isCorrectType);
            }).Repeat.Times(3).Return(true);

            map.Expect(m => m.Layers).Return(layers).Repeat.Twice();
            map.Expect(m => m.BringToFront(null)).IgnoreArguments()
                .WhenCalled(m => Assert.IsInstanceOf<IGridEditorStateRenderer>(m.Arguments[0]));

            layers.Expect(l => l.Insert(0, null)).IgnoreArguments().WhenCalled(m =>
            {
                Assert.AreEqual(0, m.Arguments[0]);
                Assert.IsInstanceOf<IGridEditorStateRenderer>(m.Arguments[1]);
            });

            layers.Expect(l => l.Remove( null)).IgnoreArguments().WhenCalled(m =>
            {
                Assert.IsInstanceOf<IGridEditorStateRenderer>(m.Arguments[0]);
            }).Return(true);

            mocks.ReplayAll();

            var gridEditorState = new GridEditorState();
            var interactor = new MapInteractor
            {
                GetCurrentMapControl = () => mapControl
            };

            // Act
            Assert.Null(interactor.PolygonMapTool);
            Assert.Null(interactor.LandBoundariesMapTool);
            Assert.Null(interactor.SplineMapTool);
            Assert.Null(interactor.Renderer);

            interactor.StartInteraction(gridEditorState);
            
            // Assert
            Assert.NotNull(interactor.PolygonMapTool);
            Assert.NotNull(interactor.LandBoundariesMapTool);
            Assert.NotNull(interactor.SplineMapTool);
            Assert.NotNull(interactor.Renderer);

            // Act
            interactor.StopInteraction();

            // Assert
            Assert.Null(interactor.PolygonMapTool);
            Assert.Null(interactor.LandBoundariesMapTool);
            Assert.Null(interactor.SplineMapTool);
            Assert.Null(interactor.Renderer);

            mocks.VerifyAll();

        }

        [Test]
        [TestCase(MapToolType.Polygon)]
        [TestCase(MapToolType.LandBoundaries)]
        [TestCase(MapToolType.None)]
        [TestCase(MapToolType.Spline)]
        [TestCase(MapToolType.InsertVertices)]
        [TestCase(MapToolType.DeleteVertices)]
        [TestCase(MapToolType.DeleteEdges)]
        [TestCase(MapToolType.MoveVertices)]
        public void GivenMapInteractor_CallingActivateTool_ShouldActiveMapToolOnMap(MapToolType toolType)
        {
            //Arrange
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var tools = mocks.StrictMock<IList<IMapTool>>();
            var map = mocks.StrictMock<IMap>();
            var layers = mocks.StrictMock<IEventedList<ILayer>>();

            mapControl.Expect(c => c.Map).Return(map).Repeat.Times(2);
            mapControl.Expect(c => c.Tools).Return(tools);
            mapControl.Expect(c => c.ActivateTool(null)).IgnoreArguments();
            mapControl.Expect(c => c.SelectTool).Return(null).Repeat.Any();
            mapControl.Expect(c => c.Refresh());

            tools.Expect(t => t.Add(null)).IgnoreArguments().Repeat.Times(7);

            map.Expect(m => m.Layers).Return(layers);
            map.Expect(m => m.BringToFront(null)).IgnoreArguments();

            layers.Expect(l => l.Insert(0, null)).IgnoreArguments();

            mocks.ReplayAll();

            // Act
            var gridEditorState = new GridEditorState();
            var interactor = new MapInteractor
            {
                GetCurrentMapControl = () => mapControl
            };
            
            interactor.StartInteraction(gridEditorState);
            interactor.ActivateTool(toolType);

            // Assert
            Assert.AreEqual(interactor.SelectedMapToolType, toolType);

            mocks.VerifyAll();
        }

        [Test]
        public void GivenMapInteractor_CallingGetUnstructuredGrids_ShouldReturnAllVisibleGridsOnTheMap()
        {
            //Arrange
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var map = mocks.StrictMock<IMap>();

            var gridLayers = new List<ILayer>
            {
                new UnstructuredGridLayer {Name = "Grid 1", Grid = new UnstructuredGrid()},
                new UnstructuredGridLayer {Name = "Grid 2", Grid = new UnstructuredGrid()},
                new VectorLayer("Other layer"),
                new UnstructuredGridLayer {Name = "Grid 3", Grid = new UnstructuredGrid()}
            };

            mapControl.Expect(c => c.Map).Return(map);
            map.Expect(m => m.GetAllVisibleLayers(false)).Return(gridLayers);

            mocks.ReplayAll();

            // Act
            var interactor = new MapInteractor
            {
                GetCurrentMapControl = () => mapControl
            };

            var grids = interactor.GetUnstructuredGrids().ToList();

            // Assert
            
            Assert.AreEqual(3, grids.Count);

            mocks.VerifyAll();
        }

        [Test]
        public void GivenMapInteractor_CallingZoomToGridExtent_ShouldZoomToMeshGeometryRendererEnvelope()
        {
            //Arrange
            var extent = new Envelope(0,10,0,20);
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var map = mocks.StrictMock<IMap>();
            var gridEditorStateRenderer = mocks.StrictMock<IGridEditorStateRenderer>();
            var disposableMeshGeometryRenderer = mocks.StrictMock<IDisposableMeshGeometryRenderer>();
            
            mapControl.Expect(c => c.Map).Return(map);
            map.Expect(m => m.ZoomToFit(extent)).IgnoreArguments();

            gridEditorStateRenderer.Expect(r => r.MeshGeometryRenderer).Return(disposableMeshGeometryRenderer);
            disposableMeshGeometryRenderer.Expect(d => d.Envelope).Return(extent);

            mocks.ReplayAll();
            
            // Act
            var interactor = new MapInteractor
            {
                GetCurrentMapControl = () => mapControl
            };
            
            TypeUtils.SetPrivatePropertyValue(interactor, nameof(interactor.Renderer), gridEditorStateRenderer);

            interactor.ZoomToGridExtent();

            // Assert
            
            mocks.VerifyAll();
        }

        [Test]
        public void GivenMapInteractor_SettingGetSplineGeometry_SetsSplineMapToolGetSplineGeometryAndRendererGetSplineGeometry()
        {
            var mocks = new MockRepository();
            var renderer = mocks.StrictMock<IGridEditorStateRenderer>();
            
            Func<Coordinate[], Coordinate[]> getSplineGeometry = (c) => c;

            renderer.Expect(r => r.SetGetSplineGeometryFunction(getSplineGeometry));
            
            mocks.ReplayAll();

            //Arrange
            var interactor = new MapInteractor();
            
            TypeUtils.SetPrivatePropertyValue(interactor, nameof(interactor.Renderer), renderer);

            // Act
            interactor.GetSplineGeometry = getSplineGeometry;

            // Assert
            mocks.VerifyAll();
        }
    }
}