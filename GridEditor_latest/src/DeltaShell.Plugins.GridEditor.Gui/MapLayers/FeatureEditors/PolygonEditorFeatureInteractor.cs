using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Drawing;
using System.Linq;
using DelftTools.Utils.Collections.Extensions;
using DelftTools.Utils.Editing;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Editors.Interactors;
using SharpMap.Styles;
using Point = NetTopologySuite.Geometries.Point;

namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.FeatureEditors
{
    public class PolygonEditorFeatureInteractor : LinearRingInteractor
    {
        private static PolygonDrawingStyle style = new PolygonDrawingStyle();

        private static readonly Bitmap TrackerSmallEnd = DrawingHelper.GenerateEllipseBitmap(style.FirstPolygonPointPen, style.EditingPointBackgroundBrush, PolygonDrawingStyle.PointSize, PolygonDrawingStyle.PointSize);
        private static readonly Bitmap TrackerSmall = DrawingHelper.GenerateEllipseBitmap(style.EditingPointBorderPen , style.EditingPointBackgroundBrush, PolygonDrawingStyle.PointSize, PolygonDrawingStyle.PointSize);
        private static readonly Bitmap SelectedTrackerSmall = DrawingHelper.GenerateEllipseBitmap(style.EditingLinePen, style.SnappingPointBackgroundBrush, PolygonDrawingStyle.PointSize, PolygonDrawingStyle.PointSize);

        public PolygonEditorFeatureInteractor(ILayer layer, IFeature feature, VectorStyle vectorStyle, IEditableObject editableObject) : base(layer, feature, vectorStyle, editableObject)
        {
            CreateTrackers();
        }

        protected override void CreateTrackers()
        {
            if (SourceFeature?.Geometry == null) return;

            Trackers.Clear();
            Trackers.AddRange(CreateTrackersForGeometry(SourceFeature.Geometry));

            AllTracker = new TrackerFeature(this, null, -1, null);
        }

        private new IEnumerable<TrackerFeature> CreateTrackersForGeometry(IGeometry geometry)
        {
            var coordinates = geometry.Coordinates; 
            if (coordinates.Length <= 2)
            {
                yield break;
            }

            for (var i = 0; i < coordinates.Length; i++)
            {
                yield return new TrackerFeature(this, new Point(coordinates[i].X, coordinates[i].Y), i, TrackerSmall);
            }

            if (coordinates.Length > 1)
            {
                yield return new TrackerFeature(this, new Point(coordinates.Last().X, coordinates.Last().Y), coordinates.Length - 1, TrackerSmallEnd);
            }
        }

        public override void SetTrackerSelection(TrackerFeature trackerFeature, bool select)
        {
            trackerFeature.Selected = select;
            trackerFeature.Bitmap = select ? SelectedTrackerSmall : TrackerSmall;
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
    }
}