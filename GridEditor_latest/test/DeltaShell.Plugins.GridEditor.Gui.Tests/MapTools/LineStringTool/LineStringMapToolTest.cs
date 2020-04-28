using System.Collections.Generic;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.UI.Forms;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool
{
    [TestFixture]
    public class LineStringMapToolTest
    {
        [Test]
        public void GivenPolygonMapTool_MouseDown_ShouldTriggerStrategyMouseDown()
        {
            //Arrange
            var envelope = new Envelope(0, 100, 0, 100);
            var worldPosition = new Coordinate(10, 10);

            var mocks = new MockRepository();
            var featureProvider = mocks.StrictMock<IFeatureProvider>();
            var layer = mocks.StrictMock<ILayer>();
            var selector = mocks.StrictMock<ISelector>();
            var map = mocks.StrictMock<IMap>();
            var mapControl = mocks.StrictMock<IMapControl>();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();

            selector.Expect(s => s.SelectedFeatureInteractors).Return(null).Repeat.AtLeastOnce();
            selector.Expect(s => s.HasSelection).Return(false).Repeat.AtLeastOnce();

            layer.Expect(l => l.GetFeatures(null)).IgnoreArguments().Return(new IFeature[0]).Repeat.Any();
            layer.Expect(l => l.Map).Return(map).Repeat.AtLeastOnce();
            layer.Expect(l => l.Envelope).Return(envelope).Repeat.AtLeastOnce();
            layer.Expect(l => l.DataSource).Return(featureProvider).Repeat.AtLeastOnce();
            
            map.Expect(m => m.PixelWidth).Return(1);
            map.Expect(m => m.UseFocusedLayer).Return(false);
            map.Expect(m => m.GetAllVisibleLayers(true)).Return(new[] {layer}).Repeat.AtLeastOnce();

            mapControl.Expect(c => c.Map).Return(map).Repeat.AtLeastOnce();

            strategy.Expect(s => s.LineStringMapToolState).Return(LineStringMapToolState.Selecting);
            strategy.Expect(s => s.HandleMouseDown(worldPosition, null, null, selector, null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                GetSelector = () => selector,
                MapControl = mapControl,
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            polygonTool.OnMouseDown(worldPosition, new MouseEventArgs(MouseButtons.Left,1, 0,0,0));

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_MapControlPreviewKeyDown_ShouldTriggerStrategyHandleMouseDown()
        {
            //Arrange
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<MapControl>();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();

            mapControl.Expect(c => c.PreviewKeyDown += null).IgnoreArguments();
            mapControl.Expect(c => c.PreviewKeyDown -= null).IgnoreArguments();

            strategy.Expect(s => s.HandleKeyDown(false, Keys.A, null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                GetSelector = () => null,
                MapControl = mapControl,
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            mapControl.Raise(c=> c.PreviewKeyDown += null, null, new PreviewKeyDownEventArgs(Keys.A));

            polygonTool.MapControl = null;

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_OnMouseDoubleClick_ShouldTriggerStrategyHandleMouseDoubleClick()
        {
            //Arrange
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();
            var layer = mocks.StrictMock<ILayer>();
            var map = mocks.StrictMock<IMap>();

            map.Expect(m => m.UseFocusedLayer).Return(false);
            map.Expect(m => m.GetAllVisibleLayers(true)).Return(new[] { layer }).Repeat.AtLeastOnce();

            mapControl.Expect(c => c.Map).Return(map).Repeat.AtLeastOnce();
            mapControl.Expect(c => c.Refresh());

            strategy.Expect(s => s.LineStringMapToolState).Return(LineStringMapToolState.Selecting);
            strategy.Expect(s => s.HandleMouseDoubleClick(null, layer, null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                GetSelector = ()=> null,
                MapControl = mapControl,
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            polygonTool.OnMouseDoubleClick(null, new MouseEventArgs(MouseButtons.Left, 0,0,0,0));

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_OnKeyUp_ShouldTriggerStrategyHandleKeyUp()
        {
            //Arrange
            var mocks = new MockRepository();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();

            strategy.Expect(s => s.HandleKeyUp(null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            polygonTool.OnKeyUp(new KeyEventArgs(Keys.K));

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_OnMouseUp_ShouldTriggerStrategyHandleMouseUp()
        {
            //Arrange
            var mocks = new MockRepository();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();
            
            strategy.Expect(s => s.LineStringMapToolState).Return(LineStringMapToolState.Selecting);
            strategy.Expect(s => s.HandleMouseUp(null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            polygonTool.OnMouseUp(new Coordinate(10,10), new MouseEventArgs(MouseButtons.Left, 0, 0, 0, 0));

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_OnMouseMove_ShouldTriggerStrategyHandleOnMouseMove()
        {
            //Arrange
            var worldPosition1 = new Coordinate(10, 10);
            var worldPosition2 = new Coordinate(20, 20);

            var mocks = new MockRepository();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();

            strategy.Expect(s => s.HandleMouseMove(null,worldPosition1,null));
            strategy.Expect(s => s.HandleMouseMove(worldPosition1, worldPosition2, null));

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                LineEditorStrategies = new List<ILineStringMapToolStrategy> { strategy }
            };

            // Act
            // first move (no previous coordinate)
            polygonTool.OnMouseMove(worldPosition1, new MouseEventArgs(MouseButtons.Left, 0, 0, 0, 0));

            // second move at same coordinate as previous coordinate (nothing happens)
            polygonTool.OnMouseMove(worldPosition1, new MouseEventArgs(MouseButtons.Left, 0, 0, 0, 0));

            // third move (real move with previous coordinate)
            polygonTool.OnMouseMove(worldPosition2, new MouseEventArgs(MouseButtons.Left, 0, 0, 0, 0));

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonMapTool_DoingMouseActionsWithOtherMouseButtons_ShouldBeIgnored()
        {
            //Arrange
            var mocks = new MockRepository();
            var strategy = mocks.StrictMock<ILineStringMapToolStrategy>();

            mocks.ReplayAll();

            var polygonTool = new LineStringMapTool(false)
            {
                LineEditorStrategies = new List<ILineStringMapToolStrategy>
                {
                    strategy
                }
            };

            // Act
            var args = new MouseEventArgs(MouseButtons.Right,0,0,0,0);

            polygonTool.OnMouseDown(new Coordinate(), args);
            polygonTool.OnMouseUp(new Coordinate(), args);
            polygonTool.OnMouseDoubleClick(null, args);

            // Assert
            mocks.VerifyAll(); // no actions should be triggered towards strategy
        }
    }
}