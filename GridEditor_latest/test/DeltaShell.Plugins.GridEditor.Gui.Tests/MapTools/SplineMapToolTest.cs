using NUnit.Framework;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools
{
    [TestFixture]
    public class SplineMapToolTest
    {
        /*[Test]
        public void GivenASetOfMouseActions_GridSplineEditorMapTool_CreatesASpline()
        {
            //Arrange
            var conversionCalled = false;
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            
            mapControl.Expect(c => c.Refresh()).Repeat.Any();
            
            mocks.ReplayAll();

            var tool = new SplineMapTool
            {
                GetSplineGeometry = coordinates =>
                {
                    conversionCalled = true;
                    return coordinates;
                },
                GridEditorState = new GridEditorState(),
                MapControl = mapControl
            };

            // Act
            using (tool)
            {
                var mouseEventArgs = new MouseEventArgs(MouseButtons.Left, 1, 0,0,0);
                tool.OnMouseDown(new Coordinate(0,0), mouseEventArgs);

                tool.OnMouseMove(new Coordinate(10, 10), mouseEventArgs);

                tool.OnMouseDown(new Coordinate(10, 10), mouseEventArgs);

                tool.OnMouseMove(new Coordinate(20, 20), mouseEventArgs);

                tool.OnMouseDoubleClick(new Coordinate(20, 20), mouseEventArgs);
            }


            // Assert

            Assert.AreEqual(1,tool.GridEditorState.Splines.Count);

            Assert.IsTrue(conversionCalled, "Coordinates should be converted towards spline coordinates");
            var expected = new []{new Coordinate(0,0), new Coordinate(10,10), new Coordinate(20,20) };
            Assert.AreEqual(expected, tool.GridEditorState.Splines[0].Geometry.Coordinates);

            mocks.VerifyAll();
        }*/
    }
}