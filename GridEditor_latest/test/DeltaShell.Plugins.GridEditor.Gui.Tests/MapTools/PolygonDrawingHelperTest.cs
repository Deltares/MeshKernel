using System;
using System.Drawing;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using DeltaShell.Plugins.GridEditor.Gui.MapTools;
using GeoAPI.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools
{
    [TestFixture]
    public class PolygonDrawingHelperTest
    {
        [Test]
        public void GivenPolygonDrawingHelper_CallingDrawLinesFromCoordinates_ShouldRenderCorrectly()
        {
            //Arrange
            var polyonCoordinates = new[]
            {
                new Coordinate(10, 10),
                new Coordinate(20, 20),
            };

            var mocks = new MockRepository();
            var map = mocks.StrictMock<IMap>();

            map.Expect(m => m.WorldToImage(null)).IgnoreArguments().Do((Func<Coordinate, PointF>) (c => new PointF((float)c.X, (float)c.Y))).Repeat.Any();
            mocks.ReplayAll();

            // Act
            using (var bitmap = new Bitmap(30,30))
            using (var graphics = Graphics.FromImage(bitmap))
            {
                PolygonDrawingHelper.DrawLinesFromCoordinates(graphics, polyonCoordinates, false, map, new PolygonDrawingStyle());

                // Assert
                var path = TestHelper.GetTestFilePath("DrawLinesFromCoordinatesReference.png");
                GuiTestHelper.CompareBitmapWithReferenceFile(bitmap, path);
            }

            mocks.VerifyAll();
        }

        [Test]
        public void GivenPolygonDrawingHelper_CallingDrawLinesFromCoordinatesWithMoreThenTwoCoordinates_ShouldRenderCorrectly()
        {
            //Arrange
            var polyonCoordinates = new[]
            {
                new Coordinate(10, 10),
                new Coordinate(20, 20),
                new Coordinate(20, 5),
            };

            var mocks = new MockRepository();
            var map = mocks.StrictMock<IMap>();

            map.Expect(m => m.WorldToImage(null)).IgnoreArguments().Do((Func<Coordinate, PointF>)(c => new PointF((float)c.X, (float)c.Y))).Repeat.Any();
            mocks.ReplayAll();

            // Act
            using (var bitmap = new Bitmap(30, 30))
            using (var graphics = Graphics.FromImage(bitmap))
            {
                PolygonDrawingHelper.DrawLinesFromCoordinates(graphics, polyonCoordinates, true, map, new PolygonDrawingStyle());

                // Assert
                var path = TestHelper.GetTestFilePath("DrawLinesFromCoordinatesThreePointsReference.png");
                GuiTestHelper.CompareBitmapWithReferenceFile(bitmap, path);
            }

            mocks.VerifyAll();
        }
    }
}