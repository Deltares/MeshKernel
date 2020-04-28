using System.Diagnostics.CodeAnalysis;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;

namespace DeltaShell.Plugins.GridEditor.Data
{
    [ExcludeFromCodeCoverage] // no logic in this class (just getters and setters)
    public class Spline : Feature
    {
        /// <summary>
        /// Geometry clicked by the user (<see cref="Feature.Geometry"/> is derived from this)
        /// </summary>
        public IGeometry UserGeometry { get; set; }

        /// <summary>
        /// Number of generated (spline) points between the <see cref="UserGeometry"/> coordinates
        /// </summary>
        public int PointsBetweenSplineCornerPoints { get; set; }
    }
}