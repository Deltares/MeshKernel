using System;
using System.Collections.Generic;
using DelftTools.TestUtils;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Extensions.CoordinateSystems;
using GeoAPI.Extensions.Coverages;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Coverages;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using SharpMap;
using SharpMap.Api;
using SharpMap.Extensions.CoordinateSystems;
using SharpMapTestUtils;

namespace DeltaShell.Plugins.GridEditor.Tests.Helpers
{
    [TestFixture]
    public class GridEditorStateImportExportExtensionsTest
    {
        [Test, Category(TestCategory.DataAccess)]
        [TestCase(ImportExportType.LandBoundary)]
        [TestCase(ImportExportType.Polygons)]
        [TestCase(ImportExportType.Spline)]
        public void GivenGridEditorState_ExportImportToShpFile_ShouldWork(ImportExportType type)
        {
            //Arrange
            var data = new GridEditorState
            {
                Polygons =
                {
                    new Feature
                    {
                        Geometry = new Polygon(new LinearRing(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(10, 0),
                            new Coordinate(0, 0),
                        }))
                    }
                },
                Splines =
                {
                    new Spline
                    {
                        UserGeometry = new LineString(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(20, 50),
                        })
                    }
                },
                LandBoundaries =
                {
                    new Feature
                    {
                        Geometry = new LineString(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(20, 50),
                        })
                    },
                    new Feature
                    {
                        Geometry = new LineString(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(20, 50),
                        })
                    }
                }

            };

            // Act

            ICoordinateSystemFactory csFactory = new OgrCoordinateSystemFactory();

            var newData = new GridEditorState{MeshCoordinateSystem = csFactory.CreateFromEPSG(28992) };
            var path = TestHelper.GetTestWorkingDirectoryGeneratedTestFilePath(".shp");

            data.ImportExport(ImportExportAction.Export, type,path);
            newData.ImportExport(ImportExportAction.Import, type, path);

            // Assert
            switch (type)
            {
                case ImportExportType.Polygons:
                    Assert.AreEqual(data.Polygons.Count, newData.Polygons.Count);
                    Assert.AreEqual(data.Polygons[0].Geometry, newData.Polygons[0].Geometry);
                    break;
                case ImportExportType.LandBoundary:
                    Assert.AreEqual(data.LandBoundaries.Count, newData.LandBoundaries.Count);
                    Assert.AreEqual(data.LandBoundaries[0].Geometry, newData.LandBoundaries[0].Geometry);
                    Assert.AreEqual(data.LandBoundaries[1].Geometry, newData.LandBoundaries[1].Geometry);
                    break;
                case ImportExportType.Spline:
                    Assert.AreEqual(data.Splines.Count, newData.Splines.Count);
                    Assert.AreEqual(data.Splines[0].UserGeometry, newData.Splines[0].UserGeometry);
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(type), type, null);
            }
        }

        [Test, Category(TestCategory.DataAccess)]
        public void GivenGridEditorState_ExportImportToNetFile_ShouldWork()
        {
            //Arrange
            var data = new GridEditorState
            {
                MeshGeometry = new DisposableMeshGeometry(UnstructuredGridTestHelper.GenerateRegularGrid(100, 100, 100, 200))
            };

            // Act
            var newData = new GridEditorState
            {
                MeshGeometry = new DisposableMeshGeometry()
            };
            var path = TestHelper.GetTestWorkingDirectoryGeneratedTestFilePath(".net");

            data.ImportExport(ImportExportAction.Export, ImportExportType.Grid, path);
            newData.ImportExport(ImportExportAction.Import, ImportExportType.Grid, path);

            // Assert
            Assert.AreEqual(data.MeshGeometry.numberOfNodes, newData.MeshGeometry.numberOfNodes);
            Assert.AreEqual(data.MeshGeometry.numberOfEdges, newData.MeshGeometry.numberOfEdges);
        }

        [Test, Category(TestCategory.DataAccess)]
        public void GivenGridEditorState_ExportImportToXYZFile_ShouldWork()
        {
            //Arrange
            var data = new GridEditorState
            {
                SamplePoints = { PointValues = new List<IPointValue>
                {
                    new PointValue{X = 10, Y = 10, value = 20},
                    new PointValue{X = 20, Y = 10, value = 40}
                }}
            };

            // Act
            var newData = new GridEditorState();
            var path = TestHelper.GetTestWorkingDirectoryGeneratedTestFilePath(".xyz");

            data.ImportExport(ImportExportAction.Export, ImportExportType.Samples, path);
            newData.ImportExport(ImportExportAction.Import, ImportExportType.Samples, path);

            // Assert
            Assert.AreEqual(data.SamplePoints.PointValues.Count, newData.SamplePoints.PointValues.Count);

            for (var index = 0; index < newData.SamplePoints.PointValues.Count; index++)
            {
                Assert.AreEqual(data.SamplePoints.PointValues[index].X, newData.SamplePoints.PointValues[index].X);
                Assert.AreEqual(data.SamplePoints.PointValues[index].Y, newData.SamplePoints.PointValues[index].Y);
                Assert.AreEqual(data.SamplePoints.PointValues[index].Value, newData.SamplePoints.PointValues[index].Value);
            }
        }

        [Test, Category(TestCategory.DataAccess)]
        public void GivenGridEditorState_ExportImportToPolFile_ShouldWork()
        {
            //Arrange
            var data = new GridEditorState
            {
                Polygons =
                {
                    new Feature
                    {
                        Geometry = new Polygon(new LinearRing(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(10, 0),
                            new Coordinate(0, 0),
                        }))
                    }
                }
            };

            // Act
            var newData = new GridEditorState();
            var path = TestHelper.GetTestWorkingDirectoryGeneratedTestFilePath(".pol");

            data.ImportExport(ImportExportAction.Export, ImportExportType.Polygons, path);
            newData.ImportExport(ImportExportAction.Import, ImportExportType.Polygons, path);

            // Assert
            Assert.AreEqual(data.Polygons.Count, newData.Polygons.Count);
            Assert.AreEqual(data.Polygons[0].Geometry, newData.Polygons[0].Geometry);
        }

        [Test, Category(TestCategory.DataAccess)]
        public void GivenGridEditorState_ExportImportToSplFile_ShouldWork()
        {
            //Arrange
            var data = new GridEditorState
            {
                Splines =
                {
                    new Spline
                    {
                        UserGeometry = new LineString(new[]
                        {
                            new Coordinate(0, 0),
                            new Coordinate(0, 10),
                            new Coordinate(10, 10),
                            new Coordinate(20, 50),
                        })
                    }
                },
            };

            // Act
            var newData = new GridEditorState();
            var path = TestHelper.GetTestWorkingDirectoryGeneratedTestFilePath(".spl");

            data.ImportExport(ImportExportAction.Export, ImportExportType.Spline, path);
            newData.ImportExport(ImportExportAction.Import, ImportExportType.Spline, path);

            // Assert
            Assert.AreEqual(data.Splines.Count, newData.Splines.Count);
            Assert.AreEqual(data.Splines[0].UserGeometry, newData.Splines[0].UserGeometry);
        }

        [Test]
        public void GivenGridEditorState_ImportingAGridWithDifferntCoordinateSystem_ShouldConvertCoordinates()
        {
            //Arrange
            var path = TestHelper.GetTestFilePath("bendprof_map.nc"); // has rd coordinates system
            ICoordinateSystemFactory csFactory = new OgrCoordinateSystemFactory();

            var state = new GridEditorState
            {
                MeshGeometry = new DisposableMeshGeometry(),
                MeshCoordinateSystem = csFactory.CreateFromEPSG(3857)
            };

            if (Map.CoordinateSystemFactory == null)
            {
                Map.CoordinateSystemFactory = new OgrCoordinateSystemFactory();
            }

            // Act
            state.ImportExport(ImportExportAction.Import, ImportExportType.Grid, path);

            // Assert
            Assert.AreEqual(368652.09377316781, state.MeshGeometry.xNodes[0]);
            Assert.AreEqual(6102279.9034865433, state.MeshGeometry.yNodes[0]);
        }
    }
}