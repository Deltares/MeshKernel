using System;
using System.Windows.Forms;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using GeoAPI.Geometries;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    public class DeleteEdgesMapTool : MapTool, IGridEditorMapTool
    {
        public GridEditorState GridEditorState { get; set; }

        public Func<Coordinate, double, bool> DeleteEdges { get; set; }

        public override void OnMouseDown(Coordinate worldPosition, MouseEventArgs e)
        {
            if (e.Button != MouseButtons.Left) return;

            // check if the coordinate to add needs to be snapped
            var coordinateToAdd = worldPosition;
            var searchSize = Math.Sqrt(Map.Envelope.Area) / 10.0;
            DeleteEdges(coordinateToAdd, searchSize);
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