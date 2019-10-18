using System.Linq;
using System.Windows.Forms;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using DeltaShell.Plugins.GridEditor.Helpers;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool.Strategies
{
    [TestFixture]
    public class AddNewPolygonLineStringMapToolStrategyTest
    {
        [Test]
        public void GivenAddNewPolygonEditorStrategy_MouseActions_ShouldResultInPolygon()
        {
            //Arrange
            var strategy = new AddNewPolygonLineStringMapToolStrategy
            {
                CommitGeometry = (g,s) => {
                    var updatedMultiPolygon = s?.GetUpdatedMultiPolygon(g as IPolygon, PolygonEditAction.Add);
                    if (updatedMultiPolygon == null) return;

                    s.Polygons.Clear();
                    s.Polygons.AddRange(updatedMultiPolygon.Geometries.Select(ge => new Feature { Geometry = ge }));
                } 
            };
            var gridState = new GridEditorState();

            // Act and Assert
            Assert.AreEqual(LineStringMapToolState.AddingNew, strategy.LineStringMapToolState);
            Assert.Null(strategy.NewGeometryCoordinates);

            // mouse down should add points
            var coordinateMouseDown = new Coordinate(0,0);
            strategy.HandleMouseDown(coordinateMouseDown, null, null, null, null );

            Assert.AreEqual(1, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseDown,strategy.NewGeometryCoordinates);

            coordinateMouseDown = new Coordinate(10, 10);
            strategy.HandleMouseDown(coordinateMouseDown, null, null, null, null);

            Assert.AreEqual(2, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseDown, strategy.NewGeometryCoordinates);

            coordinateMouseDown = new Coordinate(24, 23);
            strategy.HandleMouseDown(coordinateMouseDown, null, null, null, null);

            Assert.AreEqual(3, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseDown, strategy.NewGeometryCoordinates);

            // pressing Ctrl + z should undo adding of last point
            strategy.HandleKeyDown(true, Keys.Z, null);

            Assert.AreEqual(2, strategy.NewGeometryCoordinates.Length);
            Assert.IsFalse(strategy.NewGeometryCoordinates.Contains(coordinateMouseDown));

            // mouse move should add temporary point (to preview result)
            var coordinateMouseMove = new Coordinate(20, 20);
            strategy.HandleMouseMove(null, coordinateMouseMove, null);

            Assert.AreEqual(3, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseMove, strategy.NewGeometryCoordinates);

            var previousCoordinateMouseMove = coordinateMouseMove;
            coordinateMouseMove = new Coordinate(30, 30);
            strategy.HandleMouseMove(previousCoordinateMouseMove, coordinateMouseMove, null);

            Assert.AreEqual(3, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseMove, strategy.NewGeometryCoordinates);
            Assert.IsTrue(!strategy.NewGeometryCoordinates.Contains(previousCoordinateMouseMove));

            // add one more point to get valid polygon
            coordinateMouseDown = new Coordinate(30, 10);
            strategy.HandleMouseDown(coordinateMouseDown, null, null, null, null);

            Assert.AreEqual(3, strategy.NewGeometryCoordinates.Length);
            Assert.Contains(coordinateMouseDown, strategy.NewGeometryCoordinates);

            // double click ends editing of current polygon (and commits it to gridState)
            var coordinateMouseDoubleClick = new Coordinate(5, 5);
            strategy.HandleMouseDoubleClick(coordinateMouseDoubleClick, null, gridState);

            Assert.IsNull(strategy.NewGeometryCoordinates);
            Assert.AreEqual(1, gridState.Polygons.Count);
            Assert.AreEqual(5, gridState.Polygons[0].Geometry.Coordinates.Length);
        }
    }
}