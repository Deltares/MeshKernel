using System;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Geometries;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    public class DeleteVerticesMapTool : MapTool, IGridEditorMapTool
    {

        public GridEditorState GridEditorState { get; set; }

        public Func<int, bool> DeleteVertex { get; set; }

        public Func<Coordinate, double, int> GetVertexIndex { get; set; }

        public override void OnMouseDown(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            // check if the coordinate to add needs to be snapped
            var coordinateToAdd = worldPosition;
            var searchSize = Math.Sqrt(Map.Envelope.Area) / 100.0;
            int existingVertexIndex = GetVertexIndex(coordinateToAdd, searchSize);
            if (existingVertexIndex != -1)
            {
                DeleteVertex(existingVertexIndex);
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