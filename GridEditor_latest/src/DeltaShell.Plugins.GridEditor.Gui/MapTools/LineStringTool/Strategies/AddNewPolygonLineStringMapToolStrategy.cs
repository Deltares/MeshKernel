using System.Linq;
using DelftTools.Utils.Collections;
using DeltaShell.Plugins.GridEditor.Data;
using GeoAPI.Geometries;
using NetTopologySuite.Geometries;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal class AddNewPolygonLineStringMapToolStrategy : AddNewLineStringMapToolStrategy
    {
        protected override void CommitCoordinates(Coordinate[] lineCoordinates, GridEditorState gridEditorState)
        {
            if (lineCoordinates.Length < 3) return;

            var polygon = new Polygon(new LinearRing(lineCoordinates.Plus(lineCoordinates[0]).ToArray()));

            CommitGeometry?.Invoke(polygon, gridEditorState);
        }
    }
}

