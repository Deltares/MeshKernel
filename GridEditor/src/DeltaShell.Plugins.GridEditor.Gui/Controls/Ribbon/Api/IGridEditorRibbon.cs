using DelftTools.Shell.Gui.Forms;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api
{
    internal interface IGridEditorRibbon : IRibbonCommandHandler
    {
        /// <summary>
        /// Stops the current grid editing
        /// </summary>
        void StopGridEditing();

        /// <summary>
        /// Resets the unstructured grid list
        /// </summary>
        void ResetUnstructuredGrids();
    }
}