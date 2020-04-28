using System.Collections.Generic;
using System.Linq;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Editors.Snapping;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool
{
    [TestFixture]
    public class LineStringMapToolCoordinateSnappingExtensionsTest
    {
        [Test]
        public void GivenCoordinate_GetSnapResult_ShouldFindCorrectSnapLocation()
        {
            //Arrange
            var mocks = new MockRepository();
            var layer = mocks.StrictMock<ILayer>();

            layer.Expect(l => l.CoordinateTransformation).Return(null);

            mocks.ReplayAll();

            var coordinate = new Coordinate(1,5);
            var feature = new Feature
            {
                Geometry = new Polygon(new LinearRing(new []
                {
                    new Coordinate(0, 0),
                    new Coordinate(0, 10),
                    new Coordinate(10, 10),
                    new Coordinate(10, 0),
                    new Coordinate(0, 0)
                }))
            };

            // Act
            var snapResult = coordinate.GetSnapResult(SnapRole.FreeAtObject, feature, 2, layer);

            // Assert
            Assert.NotNull(snapResult);
            Assert.AreEqual(new Coordinate(0,5), snapResult.Location);

            mocks.VerifyAll();
        }

        [Test]
        public void GivenCoordinate_SnappedTrackerAtCoordinate_ShouldReturnCorrectTracker()
        {
            //Arrange
            var feature = new Feature
            {
                Geometry = new Polygon(new LinearRing(new[]
                {
                    new Coordinate(0, 0),
                    new Coordinate(0, 10),
                    new Coordinate(10, 10),
                    new Coordinate(10, 0),
                    new Coordinate(0, 0)
                }))
            };

            var mocks = new MockRepository();
            var map = mocks.StrictMock<IMap>();
            var layer = mocks.StrictMock<ILayer>();
            var interactor = mocks.StrictMock<IFeatureInteractor>();

            var trackers = feature.Geometry.Coordinates.Select(c => new TrackerFeature(interactor, new Point(c), 0, null)).ToList();

            interactor.Expect(i => i.Layer).Return(layer).Repeat.AtLeastOnce();
            interactor.Expect(i => i.SourceFeature).Return(feature).Repeat.AtLeastOnce();
            interactor.Expect(i => i.Trackers).Return(trackers).Repeat.AtLeastOnce();

            map.Expect(m => m.PixelWidth).Return(1);

            layer.Expect(l => l.CoordinateTransformation).Return(null);
            layer.Expect(l => l.Map).Return(map);

            mocks.ReplayAll();

            var coordinate = new Coordinate(1, 8);

            // Act
            var tracker = coordinate.SnappedTrackerAtCoordinate(interactor);

            // Assert
            Assert.NotNull(tracker);
            Assert.AreEqual(new Coordinate(0, 10), tracker.Geometry.Coordinate);

            mocks.VerifyAll();
        }

        [Test]
        public void GivenCoordinate_FindNearestFeature_ShouldReturnCorrectFeature()
        {
            //Arrange
            var feature1 = new Feature
            {
                Geometry = new Polygon(new LinearRing(new[]
                {
                    new Coordinate(0, 0),
                    new Coordinate(0, 10),
                    new Coordinate(10, 10),
                    new Coordinate(10, 0),
                    new Coordinate(0, 0)
                }))
            };

            var feature2 = new Feature
            {
                Geometry = new Polygon(new LinearRing(new[]
                {
                    new Coordinate(100, 100),
                    new Coordinate(100, 110),
                    new Coordinate(110, 110),
                    new Coordinate(110, 100),
                    new Coordinate(100, 100)
                }))
            };

            var envelope = new Envelope();
            feature1.Geometry.Coordinates.Concat(feature2.Geometry.Coordinates).ForEach(c => envelope.ExpandToInclude(c));

            var mocks = new MockRepository();
            var map = mocks.StrictMock<IMap>();
            var layer = mocks.StrictMock<ILayer>();
            var featureProvider = mocks.StrictMock<IFeatureProvider>();
            var featureRenderer = mocks.StrictMock<IFeatureRenderer>();

            map.Expect(m => m.PixelWidth).Return(1);

            featureRenderer.Expect(r => r.GetRenderedFeatureGeometry(feature1, layer)).Return(feature1.Geometry);
            featureRenderer.Expect(r => r.GetRenderedFeatureGeometry(feature2, layer)).Return(feature2.Geometry);

            layer.Expect(l => l.CustomRenderers).Return(new List<IFeatureRenderer>{ featureRenderer }).Repeat.AtLeastOnce();
            layer.Expect(l => l.DataSource).Return(featureProvider);
            layer.Expect(l => l.Envelope).Return(envelope).Repeat.AtLeastOnce();
            layer.Expect(l => l.Map).Return(map).Repeat.AtLeastOnce();
            layer.Expect(l => l.GetFeatures(new Envelope(93, 109, 100, 116))).Return(new []{feature2, feature1});

            mocks.ReplayAll();

            var coordinate = new Coordinate(101, 108);

            // Act
            var foundFeature = coordinate.FindNearestFeature(layer);

            // Assert
            Assert.AreEqual(feature2, foundFeature);

            mocks.VerifyAll();
        }
    }
}