using System.Collections.Generic;
using System.Linq;
using GeoAPI.Extensions.Feature;
using SharpMap.Api.Editors;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    /// <summary>
    /// Facade for <see cref="SelectTool"/> to hide real logic and limit coupling
    /// </summary>
    internal class SelectToolSelectorFacade : ISelector
    {
        private readonly SelectTool selectTool;

        public SelectToolSelectorFacade(SelectTool selectTool)
        {
            this.selectTool = selectTool;
        }

        /// <inheritdoc/>
        public bool HasSelection
        {
            get { return selectTool.SelectedFeatureInteractors.Any(); }
        }

        /// <inheritdoc/>
        public IEnumerable<IFeatureInteractor> SelectedFeatureInteractors
        {
            get { return selectTool.SelectedFeatureInteractors; }
        }

        /// <inheritdoc/>
        public void Select(IFeature feature)
        {
            selectTool.Select(feature);
        }

        /// <inheritdoc/>
        public void ClearSelection()
        {
            selectTool.Clear();
        }
    }
}