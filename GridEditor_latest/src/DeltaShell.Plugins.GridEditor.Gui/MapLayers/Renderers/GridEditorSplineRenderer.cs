using System.Collections.Generic;
using System.Drawing;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.CoordinateSystems.Transformations;
using SharpMap.Rendering;
using Point = NetTopologySuite.Geometries.Point;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers
{
    public class GridEditorSplineRenderer : IFeatureRenderer
    {
        static readonly IGeometryFactory GeometryFactory = new NetTopologySuite.Geometries.GeometryFactory();

        private readonly Brush splineVerticesBrush = new SolidBrush(Color.Purple);
        private readonly Pen splinePen = new Pen(Color.FromArgb(125, Color.DarkBlue), 2);

        /// <inheritdoc />
        public bool Render(IFeature feature, Graphics g, ILayer layer)
        {
            var spline = (Spline) feature;

            var lineGeometry = GetTransformedGeometry(layer, spline.Geometry) as ILineString;
            var lineUserGeometry = GetTransformedGeometry(layer, spline.UserGeometry) as ILineString;
            if (lineGeometry == null || lineUserGeometry == null) return false;

            VectorRenderingHelper.DrawLineString(g, lineGeometry, splinePen, layer.Map);

            foreach (var coordinate in lineUserGeometry.Coordinates)
            {
                VectorRenderingHelper.DrawCircle(g, new Point(coordinate.X, coordinate.Y), 5, splineVerticesBrush, layer.Map);
            }

            return true;
        }

        /// <inheritdoc />
        public IGeometry GetRenderedFeatureGeometry(IFeature feature, ILayer layer)
        {
            return GetTransformedGeometry(layer, feature.Geometry);
        }

        /// <inheritdoc />
        public IEnumerable<IFeature> GetFeatures(IGeometry geometry, ILayer layer)
        {
            return GetFeatures(geometry.EnvelopeInternal, layer);
        }

        /// <inheritdoc />
        public IEnumerable<IFeature> GetFeatures(Envelope box, ILayer layer)
        {
            return layer.GetFeatures(GeometryFactory.ToGeometry(box), false);
        }

        private static IGeometry GetTransformedGeometry(ILayer layer, IGeometry geometry)
        {
            return layer.CoordinateTransformation != null
                ? GeometryTransform.TransformGeometry(geometry, layer.CoordinateTransformation.MathTransform)
                : geometry;
        }
    }
}
