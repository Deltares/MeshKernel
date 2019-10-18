using System.Collections.Generic;
using System.Windows.Forms;
using DelftTools.Utils.Reflection;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool.Strategies
{
    [TestFixture]
    public class EditingLineStringMapToolStrategyTest
    {
        [Test]
        public void GivenEditingPolygonEditorStrategy_PressingDeleteWithSnappingResult_ShouldRemoveVertex()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var layer = mocks.StrictMock<ILayer>();

            layer.Expect(l => l.ReadOnly).Return(true);
            layer.Expect(l => l.ReadOnly = true);

            featureInteractor.Expect(i => i.Layer).Return(layer).Repeat.AtLeastOnce();
            featureInteractor.Expect(i => i.Trackers).Return(new List<TrackerFeature>
            {
                new TrackerFeature(featureInteractor,null,0,null),
                new TrackerFeature(featureInteractor, null, 1, null)
            }).Repeat.AtLeastOnce();

            featureInteractor.Expect(i => i.Start());
            featureInteractor.Expect(i => i.RemoveTracker(null)).IgnoreArguments().Return(true);
            featureInteractor.Expect(i => i.Stop());

            mocks.ReplayAll();

            var snapResult = new SnapResult(new Coordinate(10,10), null, null, null, 1,1 );

            var strategy = new EditingLineStringMapToolStrategy();

            TypeUtils.SetPrivatePropertyValue(strategy, nameof(strategy.SnapResult), snapResult);

            // Act
            strategy.HandleKeyDown(false, Keys.Delete, featureInteractor);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenEditingPolygonEditorStrategy_HandleKeyUp_ShouldRestoreReadOnlyStatusOfLayer()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var layer = mocks.StrictMock<ILayer>();

            layer.Expect(l => l.ReadOnly).Return(false);
            layer.Expect(l => l.ReadOnly = true);
            layer.Expect(l => l.ReadOnly = false);

            featureInteractor.Expect(i => i.Layer).Return(layer).Repeat.AtLeastOnce();
            featureInteractor.Expect(i => i.Trackers).Return(new List<TrackerFeature>
            {
                new TrackerFeature(featureInteractor,null,0,null),
                new TrackerFeature(featureInteractor, null, 1, null)
            }).Repeat.AtLeastOnce();

            featureInteractor.Expect(i => i.Start());
            featureInteractor.Expect(i => i.RemoveTracker(null)).IgnoreArguments().Return(true);
            featureInteractor.Expect(i => i.Stop());

            mocks.ReplayAll();

            var snapResult = new SnapResult(new Coordinate(10, 10), null, null, null, 1, 1);

            var strategy = new EditingLineStringMapToolStrategy();

            TypeUtils.SetPrivatePropertyValue(strategy, nameof(strategy.SnapResult), snapResult);

            // Act
            strategy.HandleKeyDown(false, Keys.Delete, featureInteractor);
            strategy.HandleKeyUp(featureInteractor);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenEditingPolygonEditorStrategy_HandleMouseDownWithSnapResult_ShouldTriggerInsertTracker()
        {
            //Arrange
            var expectedCoordinates = new[] { new Coordinate(10, 10), new Coordinate(20, 20) };

            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var layer = mocks.StrictMock<ILayer>();
            var feature = mocks.StrictMock<IFeature>();

            feature.Expect(f => f.Geometry).Return(new LineString(expectedCoordinates));

            layer.Expect(l => l.RenderRequired = true);

            featureInteractor.Expect(i => i.SourceFeature).Return(feature);
            featureInteractor.Expect(i => i.Layer).Return(layer).Repeat.AtLeastOnce();
            
            featureInteractor.Expect(i => i.Start());
            featureInteractor.Expect(i => i.InsertTracker(null)).IgnoreArguments().Return(true);
            featureInteractor.Expect(i => i.Stop());

            mocks.ReplayAll();

            var snapResult = new SnapResult(new Coordinate(10, 10), null, null, null, 1, 1);

            var strategy = new EditingLineStringMapToolStrategy();
            TypeUtils.SetPrivatePropertyValue(strategy, nameof(strategy.SnapResult), snapResult);
            
            // Act
            strategy.HandleMouseDown(new Coordinate(10,10),featureInteractor,null,null,null);

            // Assert
            Assert.AreEqual(expectedCoordinates, strategy.NewGeometryCoordinates);
            mocks.VerifyAll();
        }

        [Test]
        public void GivenEditingPolygonEditorStrategy_HandleMouseMove_ShouldSetSnapResultToTrackersOrPolygon()
        {
            //Arrange
            var coordinates = new[] { new Coordinate(10, 10), new Coordinate(20, 20), new Coordinate(30,10), new Coordinate(10,10) };

            var mocks = new MockRepository();
            var map = mocks.StrictMock<IMap>();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var layer = mocks.StrictMock<ILayer>();
            var feature = mocks.StrictMock<IFeature>();

            feature.Expect(f => f.Geometry).Return(new Polygon(new LinearRing(coordinates))).Repeat.AtLeastOnce();;

            map.Expect(m => m.PixelWidth).Return(1).Repeat.AtLeastOnce();

            layer.Expect(l => l.Map).Return(map);
            layer.Expect(l => l.CoordinateTransformation).Return(null);
            
            featureInteractor.Expect(i => i.SourceFeature).Return(feature).Repeat.AtLeastOnce();
            featureInteractor.Expect(i => i.Layer).Return(layer).Repeat.AtLeastOnce();

            mocks.ReplayAll();

            var snapResult = new SnapResult(new Coordinate(10, 10), null, null, null, 1, 1);

            var strategy = new EditingLineStringMapToolStrategy();
            TypeUtils.SetPrivatePropertyValue(strategy, nameof(strategy.SnapResult), snapResult);

            // Act
            strategy.HandleMouseMove(null, new Coordinate(11, 11), featureInteractor);

            Assert.AreEqual(new Coordinate(10,10), strategy.SnapResult.Location);

            // Assert
            Assert.AreEqual(coordinates, strategy.NewGeometryCoordinates);
            mocks.VerifyAll();
        }
    }
}