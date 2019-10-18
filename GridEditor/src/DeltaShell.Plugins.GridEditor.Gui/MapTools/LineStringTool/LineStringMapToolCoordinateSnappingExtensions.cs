using System;
using System.Collections.Generic;
using System.Linq;
using GeoAPI.Extensions.Feature;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;
using SharpMap.Api.Editors;
using SharpMap.Api.Layers;
using SharpMap.Editors.Snapping;
using SharpMap.Rendering;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool
{
    internal static class LineStringMapToolCoordinateSnappingExtensions
    {
        private const int PixelSnappingGravity = 16;

        public static IFeature FindNearestFeature(this Coordinate worldPos, ILayer layer)
        {
            if (layer?.DataSource == null ||
                layer.Envelope == null ||
                layer.Map == null ||
                !layer.Envelope.Contains(worldPos)) return null;

            var limit = PixelSnappingGravity * layer.Map.PixelWidth;
            var coordinateEnvelope = MapHelper.GetEnvelope(worldPos, Convert.ToSingle(limit));

            var list = layer.GetFeatures(coordinateEnvelope)?.ToList();
            if (list == null) return null;

            if (list.Count == 1)
            {
                return list[0];
            }

            var point = new Point(worldPos);
            IFeature nearest = null;
            var distance = double.MaxValue;

            foreach (var feature in list)
            {
                var geometry = layer.CustomRenderers[0].GetRenderedFeatureGeometry(feature, layer);
                var localDistance = geometry.Distance(point);

                if (localDistance >= distance || localDistance >= limit) continue;

                nearest = feature;
                distance = localDistance;
            }

            return nearest;
        }

        public static TrackerFeature SnappedTrackerAtCoordinate(this Coordinate worldPosition, IFeatureInteractor featureInteractor)
        {
            if (featureInteractor == null) return null;

            var maxDistance = PixelSnappingGravity * featureInteractor.Layer.Map.PixelWidth;
            var index = worldPosition.GetTrackerSnapResult(maxDistance, featureInteractor.SourceFeature, featureInteractor.Layer)?.SnapIndexPrevious;
            return index.HasValue ? featureInteractor.Trackers[index.Value] : null;
        }

        public static SnapResult GetTrackerSnapResult(this Coordinate worldPosition, double maxDistance, IFeature feature, ILayer layer)
        {
            var snapResultTrackers = worldPosition.GetSnapResult(SnapRole.AllTrackers, feature, maxDistance, layer);

            return snapResultTrackers != null &&
                   (snapResultTrackers.Location.Distance(worldPosition) < maxDistance)
                ? snapResultTrackers
                : null;
        }

        public static SnapResult GetSnapResult(this Coordinate worldPosition, SnapRole role, IFeature feature, double maxDistance, ILayer layer)
        {
            var snapRule = new SnapRule
            {
                Obligatory = false,
                SnapRole = role,
                PixelGravity = PixelSnappingGravity
            };

            var candidates = new[]
            {
                new Tuple<IFeature, ILayer>(feature, layer), 
            };

            var envelope = MapHelper.GetEnvelope(worldPosition, (float) maxDistance);

            return snapRule.Execute(feature, candidates, feature.Geometry, new List<IFeature> { feature },worldPosition, envelope, -1);
        }
    }
}