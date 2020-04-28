using System;
using System.Drawing;
using System.Linq;
using System.Windows;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api;
using GeoAPI.CoordinateSystems.Transformations;
using GeoAPI.Extensions.CoordinateSystems;
using GeoAPI.Geometries;
using SharpMap.Api;
using SharpMap.Api.Layers;
using SharpMap.Layers;
using SharpMap.Rendering;
using SharpMap.Rendering.Thematics;
using SharpMap.Styles;
using Point = System.Windows.Point;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers
{
    public class DisposableMeshGeometryRenderer : Layer, IDisposableMeshGeometryRenderer
    {
        private bool showPoints;
        private float pointSize = 5;
        private readonly IPrimitivesRenderer renderer = new OpenGLTaoRenderer();
        private bool drawSelectedVertices;
        private bool drawEdgeValues;
        private Envelope envelope;
        private Coordinate[] convertedVertices;
        private Func<DisposableMeshGeometry> getMeshGeometry;

        /// <inheritdoc />
        public override GeometryType GeometryType
        {
            get { return GeometryType.Unknown; }
        }

        /// <summary>
        /// Function for getting the mesh instance
        /// </summary>
        public Func<DisposableMeshGeometry> GetMeshGeometry
        {
            get { return getMeshGeometry; }
            set
            {
                getMeshGeometry = value;
                envelope = null;
            }
        }

        /// <summary>
        /// Function for getting the selected vertices
        /// </summary>
        public Func<int[]> GetSelectedVertices { get; set; }


        /// <summary>
        /// Get values on the edges (for example orthogonality value, but can also be something else)
        /// </summary>
        public Func<DisposableGeometryList> GetEdgeValues { get; set; }

        /// <summary>
        /// Function for getting coordinate system
        /// </summary>
        public Func<ICoordinateSystem> GetCoordinateSystem { get; set; }

        /// <inheritdoc />
        public bool ShowGridVertices
        {
            get { return showPoints; }
            set
            {
                showPoints = value;
                RenderRequired = true;
            }
        }

        /// <inheritdoc />
        public bool DrawEdgeValues
        {
            get { return drawEdgeValues; }
            set
            {
                drawEdgeValues = value;
                RenderRequired = true;
            }
        }

        /// <inheritdoc cref="ILayer" />
        public override Envelope Envelope
        {
            get { return envelope ?? (envelope = CalculateEnvelope()); }
        }

        /// <inheritdoc />
        public float PointSize
        {
            get { return pointSize; }
            set
            {
                pointSize = value;
                RenderRequired = true;
            }
        }

        /// <inheritdoc />
        public bool DrawSelectedVertices
        {
            get { return drawSelectedVertices; }
            set
            {
                drawSelectedVertices = value;
                RenderRequired = true;
            }
        }

        /// <inheritdoc />
        public override ICoordinateSystem CoordinateSystem
        {
            get { return GetCoordinateSystem?.Invoke(); }
        }
        
        /// <inheritdoc />
        public override ICoordinateTransformation CoordinateTransformation
        {
            get { return base.CoordinateTransformation; }
            set
            {
                base.CoordinateTransformation = value;
                MeshChanged();
            }
        }

        private Coordinate[] VertexCoordinates
        {
            get { return convertedVertices ?? (convertedVertices = GetTransformedCoordinates(GetMeshGeometry?.Invoke())); }
        }

        /// <inheritdoc />
        public override void OnRender(Graphics g, IMap map)
        {
            base.OnRender(g, map);

            var disposableMeshGeometry = GetMeshGeometry?.Invoke();
            if (disposableMeshGeometry == null || disposableMeshGeometry.numberOfNodes == 0) return;

            var worldBounds = new Rect(Map.WorldLeft,
                (Map.WorldTop - Map.WorldHeight), //bottom
                Map.Zoom, //world width
                Map.WorldHeight);

            renderer.BeginDraw(g, Map.Size, worldBounds);
            RenderEdges(disposableMeshGeometry);
            renderer.EndDraw();

            renderer.BeginDraw(g, Map.Size, worldBounds);
            RenderPoints();
            RenderSelectedVertices();
            RenderEdgeValues(disposableMeshGeometry);
            renderer.EndDraw();
        }

        /// <inheritdoc />
        public void MeshChanged()
        {
            // reset mesh related caching
            envelope = null;
            convertedVertices = null;
            RenderRequired = true;
        }

        private void RenderEdges(DisposableMeshGeometry disposableMeshGeometry)
        {
            renderer.SetDrawColor(Color.Black);

            for (int j = 0; j < disposableMeshGeometry.edgeNodes.Length; j++)
            {
                var fromIndex = disposableMeshGeometry.edgeNodes[j++];
                var toIndex = disposableMeshGeometry.edgeNodes[j];

                if (fromIndex < 0 || toIndex < 0) continue;

                var xFrom = VertexCoordinates[fromIndex].X;
                var yFrom = VertexCoordinates[fromIndex].Y;

                var xTo = VertexCoordinates[toIndex].X;
                var yTo = VertexCoordinates[toIndex].Y;

                renderer.DrawLine(xFrom, yFrom, xTo, yTo);
            }
        }

        private Coordinate[] GetTransformedCoordinates(DisposableMeshGeometry disposableMeshGeometry)
        {
            var coordinates = new Coordinate[disposableMeshGeometry.numberOfNodes];

            for (int i = 0; i < disposableMeshGeometry.numberOfNodes; i++)
            {
                var xNode = disposableMeshGeometry.xNodes[i];
                var yNode = disposableMeshGeometry.yNodes[i];

                if (xNode == -999.0 || yNode == -999.0)
                {
                    coordinates[i] = new Coordinate(-999.0, -999.0);
                    continue;
                }

                var transform = CoordinateTransformation?.MathTransform?.Transform(new []{xNode, yNode});

                coordinates[i] = transform == null
                    ? new Coordinate(xNode, yNode)
                    : new Coordinate(transform[0], transform[1]);
            }

            return coordinates;
        }

        private void RenderPoints()
        {
            if (!ShowGridVertices) return;

            renderer.SetDrawColor(Color.Red);

            for (int i = 0; i < VertexCoordinates.Length; i++)
            {
                if (VertexCoordinates[i].X == -999.0 || VertexCoordinates[i].Y == -999.0) continue;

                renderer.FillCircle(new Point(VertexCoordinates[i].X, VertexCoordinates[i].Y), PointSize);
            }
        }

        private void RenderSelectedVertices()
        {
            if (!DrawSelectedVertices) return;

            var selectedVertices = GetSelectedVertices?.Invoke();

            if (selectedVertices == null || selectedVertices.Length <= 0) return;

            renderer.SetDrawColor(Color.Green);
            
            for (int i = 0; i < selectedVertices.Length; i++)
            {
                var index = selectedVertices[i];
                var xNode = VertexCoordinates[index].X;
                var yNode = VertexCoordinates[index].Y;

                if (xNode == -999.0 || yNode == -999.0) continue;

                renderer.FillCircle(new Point(xNode, yNode), PointSize);
            }
        }
        
        private void RenderEdgeValues(DisposableMeshGeometry disposableMeshGeometry)
        {
            if (!DrawEdgeValues) return;

            var disposableGeometryListEdgeValues = GetEdgeValues?.Invoke();

            if (disposableGeometryListEdgeValues == null || disposableGeometryListEdgeValues.NumberOfCoordinates <= 0) return;

            var validEdgeZValuesArray = disposableGeometryListEdgeValues.ZCoordinates.Where(v => v > -998.999);

            if (!validEdgeZValuesArray.Any()) return;

            double minValue = validEdgeZValuesArray.Min();
            double maxValue = validEdgeZValuesArray.Max();

            var previousTheme = theme;
            var defaultStyle = new VectorStyle { GeometryType = typeof(IPoint) };
            theme = ThemeFactory.CreateGradientTheme("Edge values", defaultStyle, ColorBlend.Rainbow5,
                (float)minValue, (float)maxValue,1, 1, false, true, 12);

            int edgeIndex = 0;
            for (int j = 0; j < disposableGeometryListEdgeValues.NumberOfCoordinates; j++)
            {

                double val = disposableGeometryListEdgeValues.ZCoordinates[j];
                var fromIndex = disposableMeshGeometry.edgeNodes[edgeIndex++];
                var toIndex = disposableMeshGeometry.edgeNodes[edgeIndex++];

                if (val <= -998.999 || fromIndex < 0 || toIndex < 0) continue;

                var xFrom = VertexCoordinates[fromIndex].X;
                var yFrom = VertexCoordinates[fromIndex].Y;

                var xTo = VertexCoordinates[toIndex].X;
                var yTo = VertexCoordinates[toIndex].Y;
                
                renderer.SetDrawColor(theme.GetFillColor(val));
                renderer.DrawLine(xFrom, yFrom, xTo, yTo);
            }

            theme = previousTheme;
        }

        private Envelope CalculateEnvelope()
        {
            var calculateEnvelope = new Envelope();

            var disposableMeshGeometry = GetMeshGeometry?.Invoke();
            if (disposableMeshGeometry == null || disposableMeshGeometry.numberOfNodes == 0) return calculateEnvelope;

            for (int i = 0; i < disposableMeshGeometry.numberOfNodes; i++)
            {
                calculateEnvelope.ExpandToInclude(disposableMeshGeometry.xNodes[i], disposableMeshGeometry.yNodes[i]);
            }

            return calculateEnvelope;
        }
    }
}