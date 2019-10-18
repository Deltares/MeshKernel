using System.Collections.Generic;
using GeoAPI.Extensions.Feature;
using SharpMap.Api.Editors;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools.LineStringTool.Strategies
{
    internal interface ISelector
    {
        /// <summary>
        /// True if any features are selected
        /// </summary>
        bool HasSelection { get; }

        /// <summary>
        /// Feature inter-actors for selected features
        /// </summary>
        IEnumerable<IFeatureInteractor> SelectedFeatureInteractors { get; }

        /// <summary>
        /// Select the provided <paramref name="feature"/>
        /// </summary>
        /// <param name="feature">Feature to select</param>
        void Select(IFeature feature);

        /// <summary>
        /// Remove all selections
        /// </summary>
        void ClearSelection();
    }
}