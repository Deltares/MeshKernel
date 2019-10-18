using System.Collections.Generic;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.Api
{
    internal interface IChangePolygonMapTool
    {
        /// <summary>
        /// Get current set of <see cref="IPolygonSelection"/>
        /// </summary>
        IEnumerable<IPolygonSelection> PolygonSelections { get; }
    }
}