using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools;
using GeoAPI.Geometries;
using NUnit.Framework;
using Rhino.Mocks;
using SharpMap.Api;
using SharpMap.UI.Forms;

namespace DeltaShell.Plugins.GridEditor.Gui.Tests.MapTools
{
    [TestFixture]
    public class InsertVerticesMapToolTest
    {
        [Test]
        public void GivenASetOfMouseActions_InsertEdgesMapToolTest_CreatesNewEdges()
        {
            //Arrange
            var extent = new Envelope(0, 10, 0, 20);
            var mocks = new MockRepository();
            var mapControl = mocks.StrictMock<IMapControl>();
            var map = mocks.StrictMock<IMap>();

            mapControl.Expect(c => c.Refresh()).Repeat.Any();
            mapControl.Expect(c => c.Map.Envelope.Area).Repeat.Any();
            mapControl.Expect(c => c.Map).Return(map);
            map.Expect(m => m.Envelope.Area).Repeat.Any();

            mocks.ReplayAll();
            int vertexCounter = 0;
            int edgeCounter = 0;

            var gridEditorState = new GridEditorState
            {
                MeshGeometry = new DisposableMeshGeometry
                {
                    xNodes = new double[2],
                    yNodes = new double[2],
                    edgeNodes = new int[2]
                }
            };

            var tool = new InsertVerticesMapTool
            {
                GridEditorState = gridEditorState,

                InsertVertex = coordinates =>
                {
                    gridEditorState.MeshGeometry.xNodes[vertexCounter] = coordinates.X;
                    gridEditorState.MeshGeometry.yNodes[vertexCounter] = coordinates.Y;
                    return vertexCounter++;
                },

                InsertEdge = (firstIndex, secondIndex) =>
                {
                    gridEditorState.MeshGeometry.edgeNodes[edgeCounter] = firstIndex;
                    edgeCounter++;
                    gridEditorState.MeshGeometry.edgeNodes[edgeCounter] = secondIndex;
                    edgeCounter++;
                    return edgeCounter;
                },

                MergeVertices = () => true,

                MapControl = mapControl

            };

            tool.GetVertexIndex = (coordinate, searchRadius) => -1;

            // Act
            using (tool)
            {
                var mouseEventArgs = new MouseEventArgs(MouseButtons.Left, 1, 0, 0, 0);
                tool.OnMouseDown(new Coordinate(0, 0), mouseEventArgs);
                tool.OnMouseMove(new Coordinate(10, 10), mouseEventArgs);
                tool.OnMouseDown(new Coordinate(10, 10), mouseEventArgs);
                tool.OnMouseDoubleClick(new Coordinate(10, 10), mouseEventArgs);
            }

            // Assert
            Assert.AreEqual(2, tool.GridEditorState.MeshGeometry.xNodes.Length);
            Assert.AreEqual(2, tool.GridEditorState.MeshGeometry.yNodes.Length);
            Assert.AreEqual(2, tool.GridEditorState.MeshGeometry.edgeNodes.Length);

            Assert.AreEqual(0.0, tool.GridEditorState.MeshGeometry.xNodes[0]);
            Assert.AreEqual(0.0, tool.GridEditorState.MeshGeometry.yNodes[0]);

            Assert.AreEqual(10.0, tool.GridEditorState.MeshGeometry.xNodes[1]);
            Assert.AreEqual(10.0, tool.GridEditorState.MeshGeometry.yNodes[1]);

            mocks.VerifyAll();
        }
    }
}
