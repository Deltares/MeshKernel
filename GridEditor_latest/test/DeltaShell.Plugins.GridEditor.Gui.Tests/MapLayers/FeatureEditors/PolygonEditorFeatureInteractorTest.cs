using DelftTools.Utils.Editing;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api.Layers;
using SharpMap.Styles;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers.FeatureEditors
{
    [TestFixture]
    class PolygonEditorFeatureInteractorTest
    {
        [Test]
        public void GivenPolygonEditorFeatureInteractor_CallingStart_ShouldUpdateTrackers()
        {
            //Arrange
            var ringCoordinate = new[]
            {
                new Coordinate(0.0, 0.0, 0.0),
                new Coordinate(10.0, 0.0, 0.0),
                new Coordinate(10.0, 10.0, 0.0),
                new Coordinate(0.0, 10.0, 0.0),
                new Coordinate(0.0, 0.0, 0.0)
            };

            var feature = new Feature
            {
                Geometry = new Polygon(new LinearRing(ringCoordinate))
            };

            var vectorStyle = new VectorStyle();
            var mocks = new MockRepository();
            var layer = mocks.StrictMock<ILayer>();
            var editableObject = mocks.StrictMock<IEditableObject>();
            
            var polygonEditorFeatureInteractor = new PolygonEditorFeatureInteractor(
                layer,
                feature,
                vectorStyle,
                editableObject);

            mocks.ReplayAll();

            // Act
            polygonEditorFeatureInteractor.Start();
            polygonEditorFeatureInteractor.Stop();

            // Assert
            Assert.AreEqual(6,polygonEditorFeatureInteractor.Trackers.Count);
            
            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[0].Geometry.Coordinate.X);
            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[0].Geometry.Coordinate.Y);

            Assert.AreEqual(10.0, polygonEditorFeatureInteractor.Trackers[1].Geometry.Coordinate.X);
            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[1].Geometry.Coordinate.Y);

            Assert.AreEqual(10.0, polygonEditorFeatureInteractor.Trackers[2].Geometry.Coordinate.X);
            Assert.AreEqual(10.0, polygonEditorFeatureInteractor.Trackers[2].Geometry.Coordinate.Y);

            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[3].Geometry.Coordinate.X);
            Assert.AreEqual(10.0, polygonEditorFeatureInteractor.Trackers[3].Geometry.Coordinate.Y);

            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[4].Geometry.Coordinate.X);
            Assert.AreEqual(0.0, polygonEditorFeatureInteractor.Trackers[4].Geometry.Coordinate.Y);
            
            mocks.VerifyAll();
        }
    }
}
