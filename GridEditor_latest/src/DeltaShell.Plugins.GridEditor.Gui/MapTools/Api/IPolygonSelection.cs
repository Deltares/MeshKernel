using GeoAPI.Extensions.Feature;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.Api
{
    internal interface IPolygonSelection
    {
        /// <summary>
        /// Selected feature with the polygon geometry
        /// </summary>
        IFeature Feature { get; set; }

        /// <summary>
        /// Selected vertex indices of the polygon coordinates
        /// </summary>
        int[] SelectedIndices { get; }
    }
}