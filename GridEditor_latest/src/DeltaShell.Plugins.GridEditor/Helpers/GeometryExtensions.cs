using System;
using System.Collections.Generic;
using System.Linq;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;

namespace DeltaShell.Plugins.GridEditor.Helpers
{
    public static class GeometryExtensions
    {
        public static MultiPolygon GetClippedPolygons(this IPolygon polygon, List<IPolygon> currentPolygons)
        {
            var currentPolygon = polygon;
            foreach (var overlappingPolygon in GetOverlappingPolygons(polygon, currentPolygons))
            {
                currentPolygon = currentPolygon?.GetClippedPolygon(overlappingPolygon).FirstOrDefault();
                if (currentPolygon == null) return null;
            }

            var polygons = currentPolygons.Concat(GetPolygons(currentPolygon)).ToArray();
            return new MultiPolygon(polygons);
        }

        internal static MultiPolygon GetSubtractedOrMergedPolygons(this IPolygon polygon, PolygonEditAction polygonEditAction,
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
                        newPolygons.AddRange(overlappingPolygon.GetCombinedPolygon(polygon));
                    }

                    if (polygonEditAction == PolygonEditAction.Subtract)
                    {
                        newPolygons.AddRange(overlappingPolygon.GetSubtractedPolygon(polygon));
                    }
                }
                catch (TopologyException)
                {
                    return null;
                }
            }

            return new MultiPolygon(newPolygons.ToArray());
        }


        private static IEnumerable<IPolygon> DoPolygonOperation(Func<IGeometry> polygonOperation)
        {
            try
            {
                return GetPolygons(polygonOperation());
            }
            catch (TopologyException)
            {
                return null;
            }
        }

        public static IEnumerable<IPolygon> GetSubtractedPolygon(this IPolygon polygon, IPolygon otherPolygon)
        {
            return DoPolygonOperation(() => polygon.Difference(otherPolygon));
        }

        public static IEnumerable<IPolygon> GetCombinedPolygon(this IPolygon polygon, IPolygon otherPolygon)
        {
            return DoPolygonOperation(() => polygon.Union(otherPolygon));
        }

        public static IEnumerable<IPolygon> GetClippedPolygon(this IPolygon polygon, IPolygon otherPolygon)
        {
            return DoPolygonOperation(() => polygon.Difference(otherPolygon));
        }

        private static List<IPolygon> GetOverlappingPolygons(IPolygon polygon, List<IPolygon> currentPolygons)
        {
            return currentPolygons.Where(p => p.Intersects(polygon)).ToList();
        }

        private static IEnumerable<IPolygon> GetPolygons(IGeometry currentPolygon)
        {
            if (currentPolygon is IPolygon polygon)
            {
                return new[] { polygon };
            }

            if (currentPolygon is IMultiPolygon multiPolygon)
            {
                return multiPolygon.Geometries.OfType<IPolygon>();
            }

            return new IPolygon[0];
        }
    }
}