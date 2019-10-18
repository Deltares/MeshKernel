using System;
using System.Collections.Generic;
using System.Linq;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;

namespace DeltaShell.Plugins.GridEditor.Helpers
{
    public static class GridEditorStateExtensions
    {
        public static IMultiPolygon GetUpdatedMultiPolygon(this GridEditorState gridEditorState, IPolygon polygon, PolygonEditAction polygonEditAction)
        {
            var currentPolygons = GetCurrentPolygons(gridEditorState).ToList();

            switch (polygonEditAction)
            {
                case PolygonEditAction.Replace: return new MultiPolygon(new[] { polygon });
                case PolygonEditAction.Add: return new MultiPolygon(currentPolygons.Plus(polygon).ToArray());
                case PolygonEditAction.AddAndClip: return GetClippedPolygons(polygon, currentPolygons);
                case PolygonEditAction.Merge:
                case PolygonEditAction.Subtract:
                    return GetSubtractedOrMergedPolygons(polygon, polygonEditAction, currentPolygons);

                default: throw new ArgumentOutOfRangeException(nameof(polygonEditAction), polygonEditAction, null);
            }
        }

        private static MultiPolygon GetClippedPolygons(IPolygon polygon, List<IPolygon> currentPolygons)
        {
            IGeometry currentPolygon = polygon;
            foreach (var overlappingPolygon in GetOverlappingPolygons(polygon, currentPolygons))
            {
                try
                {
                    currentPolygon = currentPolygon.Difference(overlappingPolygon);
                }
                catch (TopologyException)
                {
                    return null;
                }
            }

            var polygons = currentPolygons.Concat(GetPolygons(currentPolygon)).ToArray();
            return new MultiPolygon(polygons);
        }

        private static MultiPolygon GetSubtractedOrMergedPolygons(IPolygon polygon, PolygonEditAction polygonEditAction,
            List<IPolygon> currentPolygons)
        {
            var overlappingPolygons = GetOverlappingPolygons(polygon, currentPolygons);
            var newPolygons = currentPolygons.Except(overlappingPolygons).ToList();

            foreach (var overlappingPolygon in overlappingPolygons)
            {
                try
                {
                    if (polygonEditAction == PolygonEditAction.Merge)
                    {
                        newPolygons.Add((IPolygon) overlappingPolygon.Union(polygon));
                    }

                    if (polygonEditAction == PolygonEditAction.Subtract)
                    {
                        newPolygons.AddRange(GetPolygons(overlappingPolygon.Difference(polygon)));
                    }
                }
                catch (TopologyException)
                {
                    return null;
                }
            }

            return new MultiPolygon(newPolygons.ToArray());
        }

        private static IEnumerable<IPolygon> GetCurrentPolygons(GridEditorState gridEditorState)
        {
            return gridEditorState?.Polygons?.Select(p => p.Geometry).OfType<IPolygon>() ?? Enumerable.Empty<IPolygon>();
        }

        private static List<IPolygon> GetOverlappingPolygons(IPolygon polygon, List<IPolygon> currentPolygons)
        {
            return currentPolygons.Where(p => p.Intersects(polygon)).ToList();
        }

        private static IEnumerable<IPolygon> GetPolygons(IGeometry currentPolygon)
        {
            if (currentPolygon is IPolygon polygon)
            {
                return new[] {polygon};
            }

            if (currentPolygon is IMultiPolygon multiPolygon)
            {
                return multiPolygon.Geometries.OfType<IPolygon>();
            }

            return new IPolygon[0];
        }
    }
}