using System;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Geometries;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    public class MoveVerticesMapTool : MapTool, IGridEditorMapTool
    {
        public GridEditorState GridEditorState { get; set; }
        public Func<Coordinate, double, int> GetVertexIndex { get; set; }
        public Func<Coordinate, int, bool> MoveVertex { get; set; }

        private int vertexIndexToMove = -1;

        public override void OnMouseDown(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            var coordinate = worldPosition;

            if (vertexIndexToMove == -1)
            {
                // check if the coordinate to add needs to be snapped
                var searchSize = Math.Sqrt(Map.Envelope.Area) / 100.0;
                vertexIndexToMove = GetVertexIndex(coordinate, searchSize);
            }
            else
            {
                MoveVertex(coordinate, vertexIndexToMove);
                vertexIndexToMove = -1;
                CommitAndRefresh();
            }
        }

        public override void OnMouseDoubleClick(object sender, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            base.OnMouseDoubleClick(sender, e);
            CommitAndRefresh();
        }

        public override void OnKeyDown(KeyEventArgs e)
        {
            if (e.KeyData == Keys.Return)
            {
                CommitAndRefresh();
                e.Handled = true;
            }
        }

        private void CommitAndRefresh()
        {
            MapControl.Refresh();
        }
    }
}