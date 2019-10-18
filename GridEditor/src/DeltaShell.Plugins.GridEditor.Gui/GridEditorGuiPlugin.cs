using System.Diagnostics;
using System.Linq;
using DelftTools.Controls;
using DelftTools.Shell.Core;
using DelftTools.Shell.Gui;
using DelftTools.Shell.Gui.Forms;
using DeltaShell.Plugins.GridEditor.Gui.Controllers;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.SharpMapGis.Gui.Forms;
using Mono.Addins;

namespace DeltaShell.Plugins.GridEditor.Gui
{
    [Extension(typeof(IPlugin))]
    public class GridEditorGuiPlugin : GuiPlugin
    {
        private IGridEditorRibbon gridEditorRibbon;
        private IGraphicsProvider graphicsProvider;

        public override string Name
        {
            get { return "GridEditorGui"; }
        }

        public override string DisplayName
        {
            get { return "Grid Editor Plugin (UI)"; }
        }

        public override string Description
        {
            get { return "Editor for creating and changing unstructured grids"; }
        }

        public override string Version
        {
            get { return GetVersionInfo().ProductVersion; }
        }

        public override string FileFormatVersion
        {
            get { return GetVersionInfo().FileVersion; }
        }

        private FileVersionInfo GetVersionInfo()
        {
            return FileVersionInfo.GetVersionInfo(GetType().Assembly.Location);
        }

        public override IGraphicsProvider GraphicsProvider
        {
            get { return graphicsProvider ?? (graphicsProvider = new GridGraphicsProvider()); }
        }

        public override IRibbonCommandHandler RibbonCommandHandler
        {
            get
            {
                return gridEditorRibbon ?? (gridEditorRibbon = new GridEditorRibbon
                {
                    MapInteractor = new MapInteractor
                    {
                        GetCurrentMapControl = () => Gui?.DocumentViews.ActiveView
                                .GetViewsOfType<MapView>()
                                .FirstOrDefault()?.MapControl
                    },
                });
            }
        }

        public override void OnActiveViewChanged(IView view)
        {
            if (!IsActive) return;

            base.OnActiveViewChanged(view);

            // ignore tool window ActiveView changes
            if (Gui.ToolWindowViews.ActiveView == view) return;

            gridEditorRibbon?.StopGridEditing();
            gridEditorRibbon?.ResetUnstructuredGrids();
        }

        public override void OnViewRemoved(IView view)
        {
            base.OnViewRemoved(view);

            var mapView = view?.GetViewsOfType<MapView>()?.FirstOrDefault();

            if (!IsActive || mapView == null) return;

            gridEditorRibbon?.StopGridEditing();
        }
    }
}