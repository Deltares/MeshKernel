using System.Collections.Generic;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool.Strategies
{
    [TestFixture]
    public class DraggingLineStringMapToolStrategyTest
    {
        [Test]
        public void GivenDraggingPolygonEditorStrategy_HandleMouseDown_ShouldStartFeatureInterActor()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var trackerFeature = new TrackerFeature(featureInteractor, new Point(0,0), 0,null);

            featureInteractor.Expect(i => i.Start());
            featureInteractor.Expect(i => i.SetTrackerSelection(trackerFeature,true));

            mocks.ReplayAll();

            var strategy = new DraggingLineStringMapToolStrategy();

            // Act
            var coordinateMouseDown = new Coordinate(0, 0);

            // mouse down should trigger start of feature inter-actor
            strategy.HandleMouseDown(coordinateMouseDown, featureInteractor, null, null, trackerFeature);

            // Assert
            mocks.VerifyAll();
        }

        [Test]
        public void GivenDraggingPolygonEditorStrategy_HandleMouseMove_ShouldTriggerFeatureInteractorMoveTracker()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();

            var trackerFeature1 = new TrackerFeature(featureInteractor, new Point(0, 0), 0, null)
            {
                Selected = true
            };

            var trackerFeature2 = new TrackerFeature(featureInteractor, new Point(3, 2), 0, null)
            {
                Selected = false
            };

            var expectedTargetFeatureCoordinates = new[]
            {
                new Coordinate(0,0),
                new Coordinate(10,10),
            };

            featureInteractor.Expect(i => i.Trackers).Return(new List<TrackerFeature> {trackerFeature1, trackerFeature2});
            featureInteractor.Expect(i => i.MoveTracker(trackerFeature1, 10, 10)).Return(true);
            featureInteractor.Expect(i => i.TargetFeature).Return(new Feature{Geometry = new LineString(expectedTargetFeatureCoordinates)});

            mocks.ReplayAll();

            var strategy = new DraggingLineStringMapToolStrategy();

            // Act
            var coordinateMouseMove = new Coordinate(10, 10);
            strategy.HandleMouseMove(new Coordinate(0,0), coordinateMouseMove, featureInteractor);

            // Assert
            Assert.AreEqual(expectedTargetFeatureCoordinates, strategy.NewGeometryCoordinates);
            mocks.VerifyAll();
        }

        [Test]
        public void GivenDraggingPolygonEditorStrategy_HandleMouseUp_ShouldTriggerFeatureInteractorStop()
        {
            //Arrange
            var mocks = new MockRepository();
            var featureInteractor = mocks.StrictMock<IFeatureInteractor>();
            var layer = mocks.StrictMock<ILayer>();

            layer.Expect(l => l.RenderRequired = true);

            featureInteractor.Expect(i => i.Stop());
            featureInteractor.Expect(i => i.Layer).Return(layer);
            
            mocks.ReplayAll();

            var strategy = new DraggingLineStringMapToolStrategy();

            // Act
            strategy.HandleMouseUp(featureInteractor);

            // Assert
            mocks.VerifyAll();
        }
    }
}