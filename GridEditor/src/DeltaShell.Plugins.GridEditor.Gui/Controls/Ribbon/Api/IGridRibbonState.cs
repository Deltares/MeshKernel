using System.Collections.Generic;
using NetTopologySuite.Extensions.Grids;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api
{
    internal interface IGridRibbonState
    {
        /// <summary>
        /// Unstructured grids that can be edited
        /// </summary>
        IList<UnstructuredGrid> UnstructuredGrids { get; set; }

        /// <summary>
        /// Selected grid to work on
        /// </summary>
        UnstructuredGrid SelectedGrid { get; set; }

        /// <summary>
        /// Refreshes the state with the current controller state
        /// </summary>
        void RefreshState();

        bool IsTaskRunning { get; set; }

        void UpdateProgress(float percentage, string progressText);
    }
}