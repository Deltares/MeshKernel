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
                case PolygonEditAction.AddAndClip: return polygon.GetClippedPolygons(currentPolygons);
                case PolygonEditAction.Merge:
                case PolygonEditAction.Subtract:
                    return polygon.GetSubtractedOrMergedPolygons(polygonEditAction, currentPolygons);

                default: throw new ArgumentOutOfRangeException(nameof(polygonEditAction), polygonEditAction, null);
            }
        }

        private static IEnumerable<IPolygon> GetCurrentPolygons(GridEditorState gridEditorState)
        {
            return gridEditorState?.Polygons?.Select(p => p.Geometry).OfType<IPolygon>() ?? Enumerable.Empty<IPolygon>();
        }
    }
}