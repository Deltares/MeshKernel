using System;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.Helpers
{
    [TestFixture]
    public class GridEditorStateExtensionsTest
    {
        [Test]
        [TestCase(PolygonEditAction.Replace)]
        [TestCase(PolygonEditAction.Add)]
        [TestCase(PolygonEditAction.AddAndClip)]
        [TestCase(PolygonEditAction.Merge)]
        [TestCase(PolygonEditAction.Subtract)]
        public void WhenNoPolygonIsPresent_CallingGetUpdatedMultiPolygonWithASpecificAction_ShouldReturnValidPolygon(PolygonEditAction action)
        {
            //Arrange
            var data = new GridEditorState();
            
            // Act
            var newPolygon = new Polygon(new LinearRing(new[]
            {
                new Coordinate(0,0),
                new Coordinate(10,0),
                new Coordinate(10,10),
                new Coordinate(0,10),
                new Coordinate(0,0),
            }));

            var newMultiPolygon = data.GetUpdatedMultiPolygon(newPolygon, action);

            // Assert
            switch (action)
            {
                case PolygonEditAction.Replace:
                    Assert.AreEqual(newPolygon, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.Add:
                    Assert.AreEqual(newPolygon, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.AddAndClip:
                    Assert.AreEqual(newPolygon, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.Merge:
                    Assert.AreEqual(0, newMultiPolygon.Geometries.Length);
                    break;
                case PolygonEditAction.Subtract:
                    Assert.AreEqual(0, newMultiPolygon.Geometries.Length);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(action), action, null);
            }
        }

        [Test]
        [TestCase(PolygonEditAction.Replace)]
        [TestCase(PolygonEditAction.Add)]
        [TestCase(PolygonEditAction.Merge)]
        [TestCase(PolygonEditAction.Subtract)]
        [TestCase(PolygonEditAction.AddAndClip)]
        public void WhenPolygonsArePresent_CallingGetUpdatedMultiPolygonWithASpecificAction_ShouldReturnValidPolygon(PolygonEditAction action)
        {
            //Arrange
            var originalPolygon = new Polygon(new LinearRing(new[]
            {
                new Coordinate(0,0),
                new Coordinate(10,0),
                new Coordinate(10,10),
                new Coordinate(0,10),
                new Coordinate(0,0),
            }));

            var data = new GridEditorState
            {
                Polygons = { new Feature { Geometry = originalPolygon } }
            };

            // Act
            var newPolygon = new Polygon(new LinearRing(new[]
            {
                new Coordinate(5,5),
                new Coordinate(15,5),
                new Coordinate(15,15),
                new Coordinate(5,15),
                new Coordinate(5,5),
            }));

            var newMultiPolygon = data.GetUpdatedMultiPolygon(newPolygon, action);

            // Assert
            switch (action)
            {
                case PolygonEditAction.Replace:
                    // replaced original polygon with new polygon
                    Assert.AreEqual(newPolygon, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.Add:
                    // added polygon as new polygon
                    Assert.AreEqual(2, newMultiPolygon.Geometries.Length);
                    Assert.AreEqual(originalPolygon, newMultiPolygon[0]);
                    Assert.AreEqual(newPolygon, newMultiPolygon[1]);
                    break;
                case PolygonEditAction.Merge:
                    // merged polygons together
                    Assert.AreEqual(1, newMultiPolygon.Geometries.Length);
                    var expectedMergedPolygons = new Polygon(new LinearRing(new[]
                    {
                        new Coordinate(10,5),
                        new Coordinate(10,0),
                        new Coordinate(0,0),
                        new Coordinate(0,10),
                        new Coordinate(5,10),
                        new Coordinate(5,15),
                        new Coordinate(15,15),
                        new Coordinate(15,5),
                        new Coordinate(10,5),
                    }));
                    Assert.AreEqual(expectedMergedPolygons, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.Subtract:
                    Assert.AreEqual(1, newMultiPolygon.Geometries.Length);
                    var expectedSubtractedPolygons = new Polygon(new LinearRing(new[]
                    {
                        new Coordinate(10,5),
                        new Coordinate(10,0),
                        new Coordinate(0,0),
                        new Coordinate(0,10),
                        new Coordinate(5,10),
                        new Coordinate(5,5),
                        new Coordinate(10,5)
                        
                    }));
                    Assert.AreEqual(expectedSubtractedPolygons, newMultiPolygon[0]);
                    break;
                case PolygonEditAction.AddAndClip:
                    Assert.AreEqual(2, newMultiPolygon.Geometries.Length);
                    Assert.AreEqual(originalPolygon, newMultiPolygon[0]);
                    var createdPolygon = new Polygon(new LinearRing(new[]
                    {
                        new Coordinate(5,10),
                        new Coordinate(5,15),
                        new Coordinate(15,15),
                        new Coordinate(15,5),
                        new Coordinate(10,5),
                        new Coordinate(10,10),
                        new Coordinate(5,10),
                    }));
                    Assert.AreEqual(createdPolygon, newMultiPolygon[1]);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(action), action, null);
            }
        }

        [Test]
        public void WhenPolygonsArePresent_CallingGetUpdatedMultiPolygonWithASubtractAction_ShouldReturnTwoPolygonIfSplittingExistingPolygon()
        {
            //            -----
            //           |*****|
            //       ----+-----+----
            //      |    |*****|    |
            //      |    |*****|    |
            //      |    |*****|    |
            //      |    |*****|    |
            //      |    |*****|    |
            //      |    |*****|    |
            //       ----+-----+----
            //           |*****|
            //            -----
            // * is new polygon

            //Arrange
            var originalPolygon = new Polygon(new LinearRing(new[]
            {
                new Coordinate(0,0),
                new Coordinate(10,0),
                new Coordinate(10,10),
                new Coordinate(0,10),
                new Coordinate(0,0),
            }));

            var data = new GridEditorState
            {
                Polygons = { new Feature { Geometry = originalPolygon } }
            };

            // Act
            var newPolygon = new Polygon(new LinearRing(new[]
            {
                new Coordinate(4,-5),
                new Coordinate(4,15),
                new Coordinate(8,15),
                new Coordinate(8,-5),
                new Coordinate(4,-5),
            }));

            var newMultiPolygon = data.GetUpdatedMultiPolygon(newPolygon, PolygonEditAction.Subtract);

            // Assert
            var expectedPolygon1 = new Polygon(new LinearRing(new []
            {
                new Coordinate(4,0),
                new Coordinate(0,0),
                new Coordinate(0,10),
                new Coordinate(4,10),
                new Coordinate(4,0),
            }));

            var expectedPolygon2 = new Polygon(new LinearRing(new[]
            {
                new Coordinate(8,10),
                new Coordinate(10,10),
                new Coordinate(10,0),
                new Coordinate(8,0),
                new Coordinate(8,10),
            }));

            Assert.AreEqual(2, newMultiPolygon.Geometries.Length);
            Assert.AreEqual(expectedPolygon1, newMultiPolygon[0]);
            Assert.AreEqual(expectedPolygon2, newMultiPolygon[1]);
        }

    }
}