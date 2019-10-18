using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using DelftTools.Utils.Collections;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Importers;
using DeltaShell.Plugins.SharpMapGis.ImportExport;
using GeoAPI.CoordinateSystems.Transformations;
using GeoAPI.Extensions.CoordinateSystems;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Extensions.Grids;
using NetTopologySuite.Geometries;
using NetTopologySuite.IO;
using SharpMap;
using SharpMap.Api;
using SharpMap.CoordinateSystems.Transformations;
using SharpMap.Data.Providers;

namespace DeltaShell.Plugins.GridEditor.Helpers
{
    public static class GridEditorStateImportExportExtensions
    {
        public static void ImportExport(this GridEditorState gridEditorState, ImportExportAction action, ImportExportType importExportType, string filePath)
        {
            var meshCoordinateSystem = gridEditorState.MeshCoordinateSystem;

            switch (action)
            {
                case ImportExportAction.Import:
                    switch (importExportType)
                    {
                        case ImportExportType.Grid:
                            ImportNetFile(gridEditorState, filePath);
                            break;
                        case ImportExportType.Polygons:

                            if (filePath.EndsWith(".shp"))
                            {
                                gridEditorState.Polygons.AddRange(GetShapeFileFeatures<Feature>(filePath, meshCoordinateSystem));
                            }
                            if (filePath.EndsWith(".pol"))
                            {
                                var createFeature = new Func<Coordinate[], Feature>(coordinates => 
                                    new Feature
                                    {
                                        Geometry = new Polygon(new LinearRing(coordinates))
                                    }
                                );
                                var polygons = InteractorFileExporterImporter.Import(filePath, createFeature);
                                gridEditorState.Polygons.AddRange(polygons);
                            }
                            break;
                        case ImportExportType.Spline:
                            if (filePath.EndsWith(".shp"))
                            {
                                var splines = GetShapeFileFeatures<Spline>(filePath, meshCoordinateSystem).ToList();
                                splines.ForEach(s =>
                                {
                                    s.UserGeometry = s.Geometry;
                                    s.Geometry = null;
                                });
                                gridEditorState.Splines.AddRange(splines);
                            }

                            if (filePath.EndsWith(".spl"))
                            {
                                var createFeature = new Func<Coordinate[], Spline>(coordinates => new Spline
                                {
                                    Geometry = null,
                                    UserGeometry = new LineString(coordinates)
                                });
                                var splines = InteractorFileExporterImporter.Import(filePath, createFeature);
                                gridEditorState.Splines.AddRange(splines);
                            }

                            break;
                        case ImportExportType.LandBoundary:
                            gridEditorState.LandBoundaries.AddRange(GetShapeFileFeatures<Feature>(filePath, meshCoordinateSystem));
                            break;
                        case ImportExportType.Samples:

                            if (filePath.EndsWith(".xyz"))
                            {
                                gridEditorState.SamplePoints.PointValues = XyzFile.Read(filePath, true); 
                            }
                            if (filePath.EndsWith(".asc"))
                            {
                                gridEditorState.SamplePoints.PointValues = AscReaderFast.Read(filePath,true);
                            }
                            break;
                        default:
                            throw new ArgumentOutOfRangeException(nameof(importExportType), importExportType, null);
                    }
                    break;
                case ImportExportAction.Export:
                    switch (importExportType)
                    {
                        case ImportExportType.Grid:
                            ExportMeshToNetFile(gridEditorState.MeshGeometry, filePath, meshCoordinateSystem);
                            break;
                        case ImportExportType.Polygons:
                            if (filePath.EndsWith(".shp"))
                            {
                                WriteShapeFileForFeatures(gridEditorState.Polygons, filePath, meshCoordinateSystem);
                            }
                            if (filePath.EndsWith(".pol"))
                            {
                                IEnumerable<IFeature> polygons = gridEditorState.Polygons;
                                InteractorFileExporterImporter.Export(ref polygons, filePath);
                            }
                            break;
                        case ImportExportType.Spline:
                            if (filePath.EndsWith(".shp"))
                            {
                                WriteShapeFileForFeatures(gridEditorState.Splines.Select(s => new Feature {Geometry = s.UserGeometry}), filePath, meshCoordinateSystem);
                            }
                            if (filePath.EndsWith(".spl"))
                            {
                                IEnumerable<IFeature> splines = gridEditorState.Splines.Select(s => new Feature { Geometry = s.UserGeometry });
                                InteractorFileExporterImporter.Export(ref splines, filePath);
                            }
                            break;
                        case ImportExportType.LandBoundary:
                            WriteShapeFileForFeatures(gridEditorState.LandBoundaries, filePath, meshCoordinateSystem);
                            break;
                        case ImportExportType.Samples:
                            XyzFile.Write(filePath, gridEditorState.SamplePoints.PointValues);
                            break;
                        default:
                            throw new ArgumentOutOfRangeException(nameof(importExportType), importExportType, null);
                    }
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(action), action, null);
            }
        }

        private static void ImportNetFile(GridEditorState gridEditorState, string filePath)
        {
            var compositeGrid = NetFileImporter.ImportGrid(filePath);

            if (compositeGrid.CoordinateSystem != null && gridEditorState.MeshCoordinateSystem != null)
            {
                var transformation = Map.CoordinateSystemFactory.CreateTransformation(compositeGrid.CoordinateSystem, gridEditorState.MeshCoordinateSystem);
                compositeGrid.Vertices.ForEach(v => TransformCoordinate(v, transformation));
            }

            compositeGrid.Vertices.AddRange(gridEditorState.MeshGeometry.CreateVertices());
            compositeGrid.Edges.AddRange(gridEditorState.MeshGeometry.CreateEdges());

            gridEditorState.MeshGeometry.Dispose();
            gridEditorState.MeshGeometry = new DisposableMeshGeometry(compositeGrid);
        }

        private static void TransformCoordinate(Coordinate coordinate, ICoordinateTransformation transformation)
        {
            coordinate.CoordinateValue = transformation.MathTransform.Transform(coordinate);
        }

        private static void ExportMeshToNetFile(DisposableMeshGeometry disposableMeshGeometry, string filePath, ICoordinateSystem meshCoordinateSystem)
        {
            var grid = new UnstructuredGrid
            {
                Vertices = disposableMeshGeometry.CreateVertices(),
                Edges = disposableMeshGeometry.CreateEdges()
            };

            NetFile.Write(filePath, grid);

            if (meshCoordinateSystem != null)
            {
                NetFile.WriteCoordinateSystem(filePath, meshCoordinateSystem);
            }
        }

        private static void WriteShapeFileForFeatures(IEnumerable<IFeature> features, string filePath, ICoordinateSystem coordinateSystem)
        {
            var geometries = features.Select(f => f.Geometry);
            ShapefileWriter.WriteGeometryCollection(filePath, new GeometryCollection(geometries.ToArray()));
            
            WriteProjFile(filePath, coordinateSystem);
        }

        private static void WriteProjFile(string path, ICoordinateSystem coordinateSystem)
        {
            if (coordinateSystem == null) return;

            var textFile = File.CreateText(path.Replace(".shp", ".prj"));
            textFile.Write(coordinateSystem.WKT);
            textFile.Close();
        }

        private static IEnumerable<T> GetShapeFileFeatures<T>(string filePath, ICoordinateSystem coordinateSystem) where T : IFeature, new()
        {
            using (var shapeFile = new ShapeFile(filePath))
            {
                var shapeFileCoordinateSystem = shapeFile.CoordinateSystem;

                return shapeFile.Features
                    .OfType<IFeature>()
                    .Select(f => new T
                    {
                        Geometry = GetTransformedGeometry(f, shapeFileCoordinateSystem, coordinateSystem),
                        Attributes = new DictionaryFeatureAttributeCollection()
                    });
            }
        }

        private static IGeometry GetTransformedGeometry(IFeature f, ICoordinateSystem shapeFileCoordinateSystem, ICoordinateSystem coordinateSystem)
        {
            if (shapeFileCoordinateSystem == null || coordinateSystem == null)
            {
                return (IGeometry) f.Geometry.Clone();
            }

            var transformation = Map.CoordinateSystemFactory.CreateTransformation(shapeFileCoordinateSystem, coordinateSystem);

            return GeometryTransform.TransformGeometry(f.Geometry, transformation.MathTransform);
        }
    }
}