using DeltaShell.Plugins.GridEditor.Data;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.Api
{
    public interface IGridEditorMapTool
    {
        /// <summary>
        /// <see cref="GridEditorState"/> to add selected items towards
        /// </summary>
        GridEditorState GridEditorState { get; set; }
    }
}