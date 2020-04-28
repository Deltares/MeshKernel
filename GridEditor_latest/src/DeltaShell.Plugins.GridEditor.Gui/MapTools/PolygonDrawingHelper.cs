using System;
using System.Drawing;
using System.Linq;
using DeltaShell.Plugins.GridEditor.Gui.MapLayers.Renderers;
using GeoAPI.Geometries;
using SharpMap.Api;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    internal static class PolygonDrawingHelper
    {
        internal static void DrawLinesFromCoordinates(Graphics graphics, Coordinate[] coordinatesToDraw, bool useLastSegmentStylePen, IMap map, PolygonDrawingStyle style)
        {
            var coordinatesToUse = coordinatesToDraw.Select(c => map.WorldToImage(c)).ToArray();

            // draw lines
            if (coordinatesToUse.Length > 1)
            {
                if (coordinatesToUse.Length > 2)
                {
                    var pointsMinusOne = new PointF[coordinatesToUse.Length - 1];
                    Array.Copy(coordinatesToUse, pointsMinusOne, coordinatesToUse.Length - 1);
                    graphics.DrawLines(style.EditingLinePen, pointsMinusOne);
                }

                var pen = useLastSegmentStylePen ? style.CreatingPolygonLastSegmentLinePen : style.EditingLinePen;
                graphics.DrawLine(pen, coordinatesToUse[coordinatesToUse.Length - 2], coordinatesToUse[coordinatesToUse.Length - 1]);
            }

            var HalfPointSize = PolygonDrawingStyle.PointSize / 2;

            // draw points 
            for (int i = 0; i < coordinatesToUse.Length; i++)
            {
                var coordinate = coordinatesToUse[i];

                var imageCoordinateX = coordinate.X - HalfPointSize;
                var imageCoordinateY = coordinate.Y - HalfPointSize;

                graphics.FillEllipse(style.EditingPointBackgroundBrush, imageCoordinateX, imageCoordinateY,
                    PolygonDrawingStyle.PointSize, PolygonDrawingStyle.PointSize);

                var pen = i == 0 ? style.FirstPolygonPointPen : style.EditingPointBorderPen;
                graphics.DrawEllipse(pen, imageCoordinateX, imageCoordinateY, PolygonDrawingStyle.PointSize,
                    PolygonDrawingStyle.PointSize);
            }
        }
    }
}