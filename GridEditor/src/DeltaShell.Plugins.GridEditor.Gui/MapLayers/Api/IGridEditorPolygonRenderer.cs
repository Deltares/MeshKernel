namespace DeltaShell.Plugins.GridEditor.Gui.MapLayers.Api
{
    public interface IGridEditorPolygonRenderer
    {
        /// <summary>
        /// Draw the polygons filled in
        /// </summary>
        bool DrawFilledInPolygons { get; set; }
    }
}