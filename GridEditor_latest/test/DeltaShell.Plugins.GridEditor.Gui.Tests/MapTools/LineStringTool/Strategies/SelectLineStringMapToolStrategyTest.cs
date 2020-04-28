using DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies;
using GeoAPI.Geometries;
using NetTopologySuite.Extensions.Features;
using NetTopologySuite.Geometries;
using NUnit.Framework;
using Rhino.Mocks;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools.LineStringTool.Strategies
{
    [TestFixture]
    public class SelectLineStringMapToolStrategyTest
    {
        [Test]
        public void GivenSelectPolygonEditorStrategy_HandleMouseDown_ShouldSelectOrDeSelectFeatures()
        {
            //Arrange
            var mocks = new MockRepository();
            var selector = mocks.StrictMock<ISelector>();

            var feature = new Feature
            {
                Geometry = new Polygon(new LinearRing(new[]
                {
                    new Coordinate(0, 0),
                    new Coordinate(0, 10),
                    new Coordinate(10, 10),
                    new Coordinate(10, 0),
                    new Coordinate(0, 0)
                }))
            };

            selector.Expect(t => t.HasSelection).Return(true);
            selector.Expect(t => t.ClearSelection());
            selector.Expect(t => t.Select(feature));

            mocks.ReplayAll();

            var strategy = new SelectLineStringMapToolStrategy();

            // Act
            var mouseDownCoordinate = new Coordinate(10, 10);
            strategy.HandleMouseDown(mouseDownCoordinate, null, feature, selector, null);

            // Assert
            Assert.AreEqual(feature.Geometry.Coordinates, strategy.NewGeometryCoordinates);

            mocks.VerifyAll();
        }
    }
}