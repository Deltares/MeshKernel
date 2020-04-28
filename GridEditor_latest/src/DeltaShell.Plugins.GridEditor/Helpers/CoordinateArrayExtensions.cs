using System.Collections.Generic;
using GeoAPI.Geometries;

namespace DeltaShell.Plugins.GridEditor.Helpers
{
    public static class CoordinateArrayExtensions
    {
        public static IEnumerable<Coordinate> RemoveDuplicateCoordinates(this IList<Coordinate> coordinateArray)
        {
            if (coordinateArray.Count > 0)
            {
                yield return coordinateArray[0];
            }

            double tolerance = 1e-3;

            for (int i = 1; i < coordinateArray.Count; i++)
            {
                var previous = coordinateArray[i - 1];
                var current = coordinateArray[i];

                if (!previous.Equals2D(current, tolerance))
                {
                    yield return current;
                }
            }
        }
    }
}