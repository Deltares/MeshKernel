using GeoAPI.Extensions.Feature;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.Api
{
    public class PolygonSelection : IPolygonSelection
    {
        /// <inheritdoc />
        public int[] SelectedIndices { get; set; }

        /// <inheritdoc />
        public IFeature Feature { get; set; }
    }
}