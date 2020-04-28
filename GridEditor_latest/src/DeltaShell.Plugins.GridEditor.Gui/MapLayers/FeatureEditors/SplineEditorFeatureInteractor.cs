using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Drawing;
using System.Linq;
using DelftTools.Utils.Collections.Extensions;
using DelftTools.Utils.Editing;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.CoordinateSystems.Transformations;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.CoordinateSystems.Transformations;
using SharpMap.Editors.Interactors;
using SharpMap.Rendering;
using SharpMap.Styles;
using Point = NetTopologySuite.Geometries.Point;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors
{
    public class SplineEditorFeatureInteractor : LinearRingInteractor
    {
        private static readonly int size = 10;
        private static readonly Bitmap TrackerVertices = DrawingHelper.GenerateEllipseBitmap(Color.Blue, Color.DarkBlue, size, size);
        private static readonly Bitmap TrackerVerticesSmall = DrawingHelper.GenerateEllipseBitmap(Color.DarkMagenta, Color.Magenta, size, size);

        private IGeometry lastSplineUserGeometry;
        
        public SplineEditorFeatureInteractor(ILayer layer, IFeature feature, VectorStyle vectorStyle,
            IEditableObject editableObject) : base(layer, feature, vectorStyle, editableObject)
        {
            CreateTrackers();
        }

        public Func<Coordinate[], Coordinate[]> GetSplineCoordinates { get; set; }

        protected override void CreateTrackers()
        {
            if (SourceFeature?.Geometry == null) return;

            Trackers.Clear();
            Trackers.AddRange(CreateTrackersForGeometry(SourceFeature));

            AllTracker = new TrackerFeature(this, null, -1, null);
        }

        private IEnumerable<TrackerFeature> CreateTrackersForGeometry(IFeature Feature)
        {
            var spline = (Spline)Feature;
            var coordinates = spline.UserGeometry?.Coordinates;

            if (coordinates == null || coordinates.Length == 0)
            {
                yield break;
            }

            for (var i = 0; i < coordinates.Length; i++)
            {
                yield return new TrackerFeature(this,
                    new Point(coordinates[i].X, coordinates[i].Y), i, TrackerVerticesSmall);
            }
        }

        public override void SetTrackerSelection(TrackerFeature trackerFeature, bool select)
        {
            trackerFeature.Selected = select;
            trackerFeature.Bitmap = select ? TrackerVerticesSmall : TrackerVertices;
        }

        public override void Stop()
        {
            var spline = SourceFeature as Spline;
            spline.UserGeometry = lastSplineUserGeometry;

            lastSplineUserGeometry = null;
            base.Stop();
        }

        public override bool MoveTracker(TrackerFeature trackerFeature, double deltaX, double deltaY, SnapResult snapResult = null)
        {
            if (trackerFeature == null)
                return false;

            var spline = SourceFeature as Spline;

            var selectedIndex = trackerFeature.Index;

            if (spline == null || selectedIndex < 0 || selectedIndex >= spline.UserGeometry.Coordinates.Length)
                return false;

            if (lastSplineUserGeometry == null)
            {
                lastSplineUserGeometry = (IGeometry)spline.UserGeometry.Clone();
            }

            var coordinate = lastSplineUserGeometry.Coordinates[selectedIndex];
            var movedCoordinate = new Coordinate(coordinate.X + deltaX, coordinate.Y + deltaY);

            trackerFeature.Geometry = new Point(movedCoordinate.Copy());
            
            lastSplineUserGeometry.Coordinates[selectedIndex].X = movedCoordinate.X;
            lastSplineUserGeometry.Coordinates[selectedIndex].Y = movedCoordinate.Y;

            TargetFeature.Geometry = new LineString(GetSplineCoordinates(lastSplineUserGeometry.Coordinates));

            return true;
        }

        public override TrackerFeature GetTrackerAtCoordinate(Coordinate worldPos)
        {
            return Trackers.FirstOrDefault(t => t.Geometry.Intersects(GetLocalExtend(t, worldPos)));
        }

        private IPolygon GetLocalExtend(TrackerFeature tracker, Coordinate worldPos)
        {
            Coordinate worldPos1 = Layer.CoordinateTransformation != null ? GeometryTransform.TransformPoint((IPoint)new Point(worldPos), this.Layer.CoordinateTransformation.MathTransform).Coordinate : worldPos;
            Coordinate coordinate = tracker.Bitmap != null ? MapHelper.ImageToWorld(this.Layer.Map, (double)tracker.Bitmap.Width, (double)tracker.Bitmap.Height) : MapHelper.ImageToWorld(this.Layer.Map, 6.0, 6.0);
            double x = coordinate.X;
            double y = coordinate.Y;
            Envelope envelope = MapHelper.GetEnvelope(worldPos1, x, y);
            IPolygon p = (IPolygon)new Polygon((ILinearRing)new LinearRing(new Coordinate[5]
            {
                new Coordinate(envelope.MinX, envelope.MinY),
                new Coordinate(envelope.MinX, envelope.MaxY),
                new Coordinate(envelope.MaxX, envelope.MaxY),
                new Coordinate(envelope.MaxX, envelope.MinY),
                new Coordinate(envelope.MinX, envelope.MinY)
            }));
            if (this.Layer.CoordinateTransformation != null)
            {
                IMathTransform transform = this.Layer.CoordinateTransformation.MathTransform.Inverse();
                p = GeometryTransform.TransformPolygon(p, transform);
            }
            return p;
        }

        [ExcludeFromCodeCoverage]
        public override bool AllowSingleClickAndMove()
        {
            return true;
        }

        [ExcludeFromCodeCoverage]
        protected override bool AllowMoveCore()
        {
            return true;
        }

        [ExcludeFromCodeCoverage]
        protected override bool AllowDeletionCore()
        {
            return true;
        }

        public override bool InsertTracker(Coordinate coordinate, SnapResult snapResult)
        {
            if (coordinate == null)
                return false;

            var spline = SourceFeature as Spline;

            if (spline == null)
                return false;

            if (lastSplineUserGeometry == null)
            {
                lastSplineUserGeometry = (IGeometry)spline.UserGeometry.Clone();
            }

            int index = snapResult.SnapIndexPrevious / spline.PointsBetweenSplineCornerPoints + 1;
            var coordinates = lastSplineUserGeometry.Coordinates.ToList();
            coordinates.Insert(index, new Coordinate(coordinate.X, coordinate.Y));
            lastSplineUserGeometry = new LineString(coordinates.ToArray());
            TargetFeature.Geometry = new LineString(GetSplineCoordinates(lastSplineUserGeometry.Coordinates));

            return true;
        }

    }
}
