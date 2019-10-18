using System;
using System.Collections.Generic;
using System.Drawing;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.CoordinateSystems.Transformations;
using SharpMap.Rendering;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers
{
    public class GridEditorPolygonRenderer : IFeatureRenderer, IGridEditorPolygonRenderer, IDisposable
    {
        static IGeometryFactory geometryFactory = new NetTopologySuite.Geometries.GeometryFactory();

        readonly PolygonDrawingStyle style = new PolygonDrawingStyle();

        public bool DrawFilledInPolygons { get; set; } = true;

        public bool Render(IFeature feature, Graphics g, ILayer layer)
        {
            var currentGeometry = layer.CoordinateTransformation != null
                ? GeometryTransform.TransformGeometry(feature.Geometry, layer.CoordinateTransformation.MathTransform)
                : feature.Geometry;

            var polygon = currentGeometry as IPolygon;

            if (polygon == null) return false;

            var brush = DrawFilledInPolygons ? style.PolygonBrush : null;

            VectorRenderingHelper.DrawPolygon(g, polygon, brush, style.PolygonBorderPen, false, layer.Map);

            for (int j = 0; j < polygon.Coordinates.Length; j++)
            {
                var coordinate = polygon.Coordinates[j];
                VectorRenderingHelper.DrawCircle(g, new NetTopologySuite.Geometries.Point(coordinate.X, coordinate.Y), 5, style.PolygonBrush, layer.Map);
            }
        
            return true;
        }

        public IGeometry GetRenderedFeatureGeometry(IFeature feature, ILayer layer)
        {
            return layer.CoordinateTransformation != null
                ? GeometryTransform.TransformGeometry(feature.Geometry, layer.CoordinateTransformation.MathTransform)
                : feature.Geometry;
        }

        public IEnumerable<IFeature> GetFeatures(IGeometry geometry, ILayer layer)
        {
            return GetFeatures(geometry.EnvelopeInternal, layer);
        }

        public IEnumerable<IFeature> GetFeatures(Envelope box, ILayer layer)
        {
            return layer.GetFeatures(geometryFactory.ToGeometry(box), false);
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            style?.Dispose();
        }
    }
}