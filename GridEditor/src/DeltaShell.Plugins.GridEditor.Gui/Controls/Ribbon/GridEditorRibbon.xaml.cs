using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using DelftTools.Controls;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.GridEditor.Helpers;
using Microsoft.Win32;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon
{
    /// <summary>
    /// Interaction logic for GridEditorRibbon.xaml
    /// </summary>
    internal partial class GridEditorRibbon : UserControl, IGridEditorRibbon
    {
        public GridEditorRibbon()
        {
            InitializeComponent();

            GridEditorTab.Group = geospatialContextualGroup;
            GridEditingGroupViewModel.Controller.CommitState = CommitState;
            GridEditingGroupViewModel.GetFilePath = GetFilePath;
        }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public IEnumerable<ICommand> Commands
        {
            get { yield break; }
        }

        /// <summary>
        /// <see cref="IMapInteractor"/> for doing map related actions
        /// </summary>
        public IMapInteractor MapInteractor
        {
            get { return GridEditingGroupViewModel.Controller.MapInteractor; }
            set { GridEditingGroupViewModel.Controller.MapInteractor = value; }
        }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public object GetRibbonControl()
        {
            return RibbonControl;
        }

        /// <inheritdoc />
        [ExcludeFromCodeCoverage]
        public void ValidateItems()
        {
            GridEditingGroupViewModel.RefreshState();
        }

        /// <inheritdoc />
        public bool IsContextualTabVisible(string tabGroupName, string tabName)
        {
            if (tabGroupName != "geospatialContextualGroup" || tabName != "GridEditorTab")
                return false;

            var unstructuredGrids = MapInteractor?.GetUnstructuredGrids();

            return unstructuredGrids != null && unstructuredGrids.Any() || GridEditingGroupViewModel.IsEditing;
        }

        /// <inheritdoc />
        public void StopGridEditing()
        {
            if (!GridEditingGroupViewModel.Controller.IsEditing) return;

            GridEditingGroupViewModel.Controller.IsEditing = false;
        }

        /// <inheritdoc />
        public void ResetUnstructuredGrids()
        {
            GridEditingGroupViewModel?.Controller?.ResetUnstructuredGrids();
        }

        [ExcludeFromCodeCoverage] // exclude because tests for this will make the build-server hang
        private bool CommitState()
        {
            var message = "Do you want to save the changes that were made to the grid";
            return MessageBox.Show(message, "Save changes", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.Yes) == MessageBoxResult.Yes;
        }

        [ExcludeFromCodeCoverage] // exclude because tests for this will make the build-server hang
        private string GetFilePath(string fileFilter, ImportExportAction action)
        {
            FileDialog dialog;
            switch (action)
            {
                case ImportExportAction.Import:
                    dialog = new OpenFileDialog { Filter = fileFilter };
                    break;
                case ImportExportAction.Export:
                    dialog = new SaveFileDialog { Filter = fileFilter };
                    break;
                default:
                    throw new ArgumentOutOfRangeException(nameof(action), action, null);
            }

            var showDialogResult = dialog.ShowDialog();

            return showDialogResult.HasValue && showDialogResult.Value
                ? dialog.FileName
                : null;
        }
    }
}
