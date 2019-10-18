using System.Collections.Generic;
using System.Linq;
using DeltaShell.Plugins.GridEditor.Gui.MapTools.Api;
using SharpMap.Api.Editors;
using SharpMap.UI.Tools;

namespace DeltaShell.Plugins.GridEditor.Gui.MapTools
{
    internal class ChangePolygonMapTool : IChangePolygonMapTool
    {
        public MoveTool MoveTool { get; } = new MoveTool
        {
            Name = "ChangePolygonMapToolMoveTool",
            FallOffPolicy = FallOffType.Ring
        };

        public IEnumerable<IPolygonSelection> PolygonSelections
        {
            get
            {
                var featureInteractors = MoveTool.MapControl?.SelectTool?.SelectedFeatureInteractors;
                if (featureInteractors == null) yield break;

                foreach (var featureInteractor in featureInteractors)
                {
                    yield return new PolygonSelection
                    {
                        Feature = featureInteractor.SourceFeature,
                        SelectedIndices = featureInteractor.Trackers
                            .Where(t => t.Selected)
                            .Select(t => t.Index)
                            .ToArray()
                    };
                }
            }
        }
    }
}