using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using SharpMap.CoordinateSystems.Transformations;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    public class InsertVerticesMapTool : MapTool, IGridEditorMapTool, IDisposable
    {
        private readonly IList<Coordinate> coordinates = new List<Coordinate>();
        private readonly IList<int> verticesIndexes = new List<int>();
        private Coordinate previousCoordinate;
        private readonly Brush edgePointBrush = new SolidBrush(Color.Red);
        private readonly Pen edgePen = new Pen(Color.Red, 2);
        private bool drawSnappingVertex = false;

        public GridEditorState GridEditorState { get; set; }
        public Func<Coordinate, int> InsertVertex { get; set; }
        public Func<int, int, int> InsertEdge { get; set; }
        public Func<bool> MergeVertices { get; set; }
        public Func<Coordinate, double, int> GetVertexIndex { get; set; }

        public override bool IsActive
        {
            get { return base.IsActive; }
            set
            {
                base.IsActive = value;
                CancelCurrentLine();
            }
        }

        protected void CancelCurrentLine()
        {
            coordinates.Clear();
            verticesIndexes.Clear();
            previousCoordinate = null;
        }

        public override void OnMouseDown(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            // check if the coordinate to add needs to be snapped
            var coordinateToAdd = worldPosition;
            var searchSize = Math.Sqrt(Map.Envelope.Area) / 100.0;
            int existingVertexIndex = GetVertexIndex(coordinateToAdd, searchSize);
            if (existingVertexIndex != -1)
            {
                coordinateToAdd.X = GridEditorState.MeshGeometry.xNodes[existingVertexIndex];
                coordinateToAdd.Y = GridEditorState.MeshGeometry.yNodes[existingVertexIndex];
            }

            coordinates.Add(coordinateToAdd);

            if (coordinates == null || coordinates.Count <= 1)
            {
                return;
            }

            var coordinateArrayFiltered = coordinates.RemoveDuplicateCoordinates()
                .ToArray();

            if (GridEditorState.MeshCoordinateSystem != null && Map.CoordinateSystem != null)
            {
                var transformation = SharpMap.Map.CoordinateSystemFactory.CreateTransformation(Map.CoordinateSystem, GridEditorState.MeshCoordinateSystem);
                coordinateArrayFiltered = GeometryTransform.TransformLineString(new LineString(coordinateArrayFiltered), transformation.MathTransform).Coordinates;
            }

            CommitCoordinates(coordinateArrayFiltered);

        }

        private void CommitCoordinates(Coordinate[] toArray)
        {
            if (toArray == null || toArray.Length < 2) return;
            if (toArray.Length == 2)
            {
                verticesIndexes.Add(InsertVertex(toArray[0]));
                verticesIndexes.Add(InsertVertex(toArray[1]));
            }
            else
            {
                verticesIndexes.Add(InsertVertex(toArray.Last()));
            }

            InsertEdge(verticesIndexes[verticesIndexes.Count - 2], verticesIndexes[verticesIndexes.Count - 1]);
        }

        public override void OnMouseMove(Coordinate worldPosition, MouseEventArgs e)
        {
            base.OnMouseMove(worldPosition, e);

            if (coordinates.Count == 0 || (previousCoordinate != null && worldPosition.Equals2D(previousCoordinate)))
            {
                return;
            }

            previousCoordinate = worldPosition;

            MapControl.Refresh();
        }

        public override void OnMouseDoubleClick(object sender, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            base.OnMouseDoubleClick(sender, e);
            MergeVertices();
            CommitAndRefresh();
        }

        public override void OnKeyDown(KeyEventArgs e)
        {
            if (e.KeyData == Keys.Return)
            {
                CommitAndRefresh();
                e.Handled = true;
            }

            if (e.KeyData == Keys.F1)
            {
                drawSnappingVertex = true;
            }
        }

        public override void OnPaint(Graphics graphics)
        {
            base.OnPaint(graphics);

            var coordinatesToUse = (previousCoordinate != null
                ? coordinates.Plus(previousCoordinate)
                : coordinates).ToArray();

            coordinatesToUse = coordinatesToUse.RemoveDuplicateCoordinates().ToArray();

            var renderingCoordinates = coordinatesToUse
                .Select(c => Map.WorldToImage(c))
                .ToArray();

            // draw line
            if (renderingCoordinates.Length >= 2)
            {
                for (int i = 0; i < renderingCoordinates.Length - 1; i++)
                {
                    graphics.DrawLine(edgePen, renderingCoordinates[i].X, renderingCoordinates[i].Y, renderingCoordinates[i + 1].X, renderingCoordinates[i + 1].Y);
                }
            }

            // draw points
            var size = 6;
            var halfSize = size / 2;
            foreach (var point in renderingCoordinates)
            {
                graphics.FillEllipse(edgePointBrush, (float)(int)(point.X - halfSize), (int)(point.Y - halfSize), size, size);
            }

            // Draw the snapping vertex if required 
            if (coordinatesToUse.Length <= 0) return;
            if (!drawSnappingVertex) return;
            var lastCoordinatesToUse = coordinatesToUse.Last();
            var lastRenderingCoordinate = renderingCoordinates.Last();
            var searchSize = Math.Sqrt(Map.Envelope.Area) / 100.0;
            if (GetVertexIndex(lastCoordinatesToUse, searchSize) != -1)
            {
                graphics.FillEllipse(edgePointBrush, (float)(int)(lastRenderingCoordinate.X - size), (int)(lastRenderingCoordinate.Y - size), size * 2, size * 2);
            }
            drawSnappingVertex = false;
        }

        private void CommitAndRefresh()
        {

            CancelCurrentLine();
            MapControl.Refresh();
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (!disposing) return;
            edgePointBrush?.Dispose();
            edgePen?.Dispose();
        }
    }
}