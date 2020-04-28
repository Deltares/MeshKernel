using System.Drawing;
using DelftTools.Utils.Editing;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Styles;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapLayers.FeatureEditors
{
    [TestFixture]
    class SplineEditorFeatureInteractorTest
    {
        [Test]
        public void GivenSplineEditorFeatureInteractor_CallingMoveTracker_ShouldUpdateTargetFeatureGeometry()
        {
            //Arrange
            var mocks = new MockRepository();

            var layer = mocks.StrictMock<ILayer>();
            var cornerPoints = new[] { new Coordinate(10.0, 10.0, 0.0), new Coordinate(50.0, 50.0, 0.0) };
            var inputLine = new LineString(cornerPoints);
            var spline = new Spline { UserGeometry = inputLine };
            var feature = spline;

            var vectorStyle = new VectorStyle();
            var editableObject = mocks.StrictMock<IEditableObject>();
            var bitmap = new Bitmap(10,10);

            var splineFeatureInteractor = new SplineEditorFeatureInteractor(layer, feature, vectorStyle, editableObject)
            {
                GetSplineCoordinates = (coord) => coord
            };

            var trackerFeature = new TrackerFeature(splineFeatureInteractor, inputLine, 0, bitmap);

            mocks.ReplayAll();
            // Act
            splineFeatureInteractor.Start();
            splineFeatureInteractor.MoveTracker(trackerFeature, 10.0, 10.0, null);
            splineFeatureInteractor.Stop();

            // Assert
            Assert.AreEqual(20.0,splineFeatureInteractor.TargetFeature.Geometry.Coordinates[0].X);
            Assert.AreEqual(20.0, splineFeatureInteractor.TargetFeature.Geometry.Coordinates[0].Y);
            Assert.AreEqual(50.0, splineFeatureInteractor.TargetFeature.Geometry.Coordinates[1].X);
            Assert.AreEqual(50.0, splineFeatureInteractor.TargetFeature.Geometry.Coordinates[1].Y);
            mocks.VerifyAll();
        }
    }
}
