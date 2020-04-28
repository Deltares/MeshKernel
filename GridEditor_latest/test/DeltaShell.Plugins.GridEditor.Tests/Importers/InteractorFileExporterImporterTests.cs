using System;
using System.Collections.Generic;
using System.IO;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Importers;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Tests.Importers
{
    [TestFixture]
    public class InteractorFileExporterImporterTests
    {
        [Test]
        public void GivenSplFile_Import_ShouldReturnCorrectCornerPoints()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("ValidSplineFile.spl");

            // Act
            var createFeature = new Func<Coordinate[], Spline>(coordinates => new Spline
            {
                Geometry = null,
                UserGeometry = new LineString(coordinates)
            });

            var splines = (List<Spline>) InteractorFileExporterImporter.Import(path, createFeature);

            // Assert
            Assert.NotNull(splines);
            Assert.AreEqual(2, splines.Count);
        }

        [Test]
        public void GivenSplFile_Export_ShouldWriteCornerPoints()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("ExportedSplineFile.spl");
            var validExportedPath = TestHelper.GetTestFilePath("ValidExportedSplineFile.spl");

            var firstSplineCoordinates = new[]
            {
                new Coordinate(10.0,10.0,0.0),
                new Coordinate(20.0,20.0,0.0)
            };

            var secondSplineCoordinates = new[]
            {
                new Coordinate(30.0,30.0,0.0),
                new Coordinate(40.0,40.0,0.0)
            };

            Spline firstSpline = new Spline
            {
                Geometry = new LineString(firstSplineCoordinates),
                UserGeometry = null
            };

            Spline secondSpline = new Spline
            {
                Geometry = new LineString(secondSplineCoordinates),
                UserGeometry = null
            };

            IEnumerable<IFeature> splines = new List<Spline>
            {
                firstSpline,
                secondSpline
            };

            // Act
            InteractorFileExporterImporter.Export(ref splines, path);
            //Assert
            using (var test = new StreamReader(path))
            using (var valid = new StreamReader(validExportedPath))
            {
                FileAssert.AreEqual(test.BaseStream, valid.BaseStream);
            }
            File.Delete(path);
        }

        //Polygon files in interactor format

        [Test]
        public void GivenPolFile_Import_ShouldWriteCornerPoints()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("ValidPolygonFile.pol");

            // Act
            var createFeature = new Func<Coordinate[], Feature>(coordinates =>
                new Feature
                {
                    Geometry = new Polygon(new LinearRing(coordinates))
                }
            );

            var polygons = (List<Feature>)InteractorFileExporterImporter.Import(path, createFeature);

            // Assert
            Assert.NotNull(polygons);
            Assert.AreEqual(2, polygons.Count);
        }

        [Test]
        public void GivenPolFile_Export_ShouldWriteCornerPoints()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("ExportedPolygonFile.pol");
            var validExportedPath = TestHelper.GetTestFilePath("ValidExportedPolygonFile.pol");

            var firstPolygonRingCoordinates = new[]
            {
                new Coordinate(393.263481920117, 1574.58227864247, 0.0),
                new Coordinate(-202.288952877188, 956.661677280008, 0.0),
                new Coordinate(829.442729940678, 584.790908134274, 0.0),
                new Coordinate(1318.74637355349, 1549.41809125667, 0.0),
                new Coordinate(393.263481920117, 1574.58227864247, 0.0)
            };

            var firstPolygonFeature = new Feature
            {
                Geometry = new Polygon(new LinearRing(firstPolygonRingCoordinates))
            };

            var secondPolygonRingCoordinates = new[]
            {
                new Coordinate(1905.91074588886, 1510.27379976764, 0.0),
                new Coordinate(885.363146353571, 143.019618472423, 0.0),
                new Coordinate(2255.41334846944, 224.104222271117, 0.0),
                new Coordinate(2579.75176366421, 1294.98019657801, 0.0),
                new Coordinate(1905.91074588886, 1510.27379976764, 0.0)
            };

            var secondPolygonFeature = new Feature
            {
                Geometry = new Polygon(new LinearRing(secondPolygonRingCoordinates))
            };

            IEnumerable<IFeature> polygons = new List<Feature>
            {
                firstPolygonFeature,
                secondPolygonFeature
            };

            // Act
            InteractorFileExporterImporter.Export(ref polygons, path);
            //Assert
            using(var test = new StreamReader(path))
            using (var valid = new StreamReader(validExportedPath))
            {
                FileAssert.AreEqual(test.BaseStream, valid.BaseStream);
            }
            File.Delete(path);
        }
    }
}
