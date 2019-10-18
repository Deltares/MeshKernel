using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Windows.Input;
using DelftTools.Controls.Wpf.Commands;
using DelftTools.Utils;
using DeltaShell.Plugins.GridEditor.Data;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controllers;
using DeltaShell.Plugins.GridEditor.Gui.Controllers.Api;
using DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.Api;
using DeltaShell.Plugins.GridEditor.Helpers;
using NetTopologySuite.Extensions.Grids;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon
{
    /// <summary>
    /// View model for <see cref="T:DeltaShell.Plugins.GridEditor.Gui.Controls.Ribbon.GridEditorRibbon" />
    /// </summary>
    /// <inheritdoc cref="IGridRibbonState"/>
    /// <inheritdoc cref="INotifyPropertyChanged"/>
    public class GridEditorRibbonViewModel : INotifyPropertyChanged, IGridRibbonState
    {
        private IList<UnstructuredGrid> unstructuredGrids;
        private UnstructuredGrid selectedGrid;
        private ImportExportType selectedImportType = ImportExportType.Grid;
        private ImportExportType selectedExportType;
        private float currentProgress;
        private string currentProgressText;
        private bool isTaskRunning;

        public GridEditorRibbonViewModel() : this(null)
        {
        }

        internal GridEditorRibbonViewModel(IGridEditorController gridEditorController)
        {
            Controller = gridEditorController ?? new GridEditorController {RibbonState = this};

            SetCommands();
        }

        /// <summary>
        /// Is in edit state (mode)
        /// </summary>
        public bool IsEditing
        {
            get { return Controller.IsEditing; }
            set
            {
                if (IsEditing == value) return;

                Controller.IsEditing = value;
                OnPropertyChanged();
                OnPropertyChanged(nameof(CanDoActions));
            }
        }

        public bool CanDoActions
        {
            get { return IsEditing && !IsTaskRunning; }
        }

        /// <summary>
        /// Selected export type
        /// </summary>
        public ImportExportType SelectedExportType
        {
            get { return selectedExportType; }
            set
            {
                selectedExportType = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Selected import type
        /// </summary>
        public ImportExportType SelectedImportType
        {
            get { return selectedImportType; }
            set
            {
                selectedImportType = value;
                OnPropertyChanged();
                OnPropertyChanged(nameof(SelectedImportTypeString));
            }
        }

        /// <summary>
        /// <see cref="SelectedImportType"/> as string (needed because of ribbon logic (two line label control))
        /// </summary>
        public string SelectedImportTypeString
        {
            get { return "Import " + SelectedImportType; }
        }

        /// <summary>
        /// <see cref="SelectedImportType"/> as string (needed because of ribbon logic (two line label control))
        /// </summary>
        public string SelectedExportTypeString
        {
            get { return "Export " + SelectedExportType; }
        }

        /// <summary>
        /// Check if there is more than one grid to select from
        /// </summary>
        public bool MoreThanOneGrid
        {
            get { return UnstructuredGrids?.Count > 1; }
        }

        /// <summary>
        /// Currently selected work flow
        /// </summary>
        public GenerateGridWorkFlowType SelectedGenerateGridWorkFlowType
        {
            get { return Controller.SelectedGenerateGridWorkFlowType; }
            set
            {
                Controller.SelectedGenerateGridWorkFlowType = value;

                OnPropertyChanged();

                // force update of related properties
                OnPropertyChanged(nameof(SelectedGenerateGridWorkFlowTypeString));
                OnPropertyChanged(nameof(IsPolygonToolVisible));
                OnPropertyChanged(nameof(IsSplineToolVisible));
            }
        }

        /// <summary>
        /// Name of the <see cref="SelectedGenerateGridWorkFlowType"/>
        /// </summary>
        public string SelectedGenerateGridWorkFlowTypeString
        {
            get { return SelectedGenerateGridWorkFlowType.ToString(); }
        }

        /// <summary>
        /// Can edit the selected grid
        /// </summary>
        public bool CanEdit
        {
            get { return SelectedGrid != null && !IsTaskRunning; }
        }

        /// <inheritdoc />
        public IList<UnstructuredGrid> UnstructuredGrids
        {
            get { return unstructuredGrids; }
            set
            {
                unstructuredGrids = value;
                if (SelectedGrid == null || unstructuredGrids != null && unstructuredGrids.Any(g => g != SelectedGrid))
                {
                    SelectedGrid = unstructuredGrids?.FirstOrDefault();
                }

                OnPropertyChanged();
                OnPropertyChanged(nameof(MoreThanOneGrid));
            }
        }

        /// <inheritdoc/>
        public UnstructuredGrid SelectedGrid
        {
            get { return selectedGrid; }
            set
            {
                selectedGrid = value;
                OnPropertyChanged();
                OnPropertyChanged(nameof(CanEdit));
            }
        }

        /// <summary>
        /// Show grid vertices (when in edit mode <see cref="IsEditing"/>)
        /// </summary>
        public bool ShowVertices
        {
            get { return Controller.ShowVertices; }
            set
            {
                Controller.ShowVertices = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Draw the <see cref="GridEditorState.Polygons"/> filled in
        /// </summary>
        public bool DrawFilledInPolygons
        {
            get { return Controller.DrawFilledInPolygons; }
            set
            {
                Controller.DrawFilledInPolygons = value;
                OnPropertyChanged();
            }
        }



        /// <summary>
        /// Draws the <see cref="GridEditorState.SelectedVertices"/>
        /// </summary>
        public bool DrawSelectedVertices
        {
            get { return Controller.DrawSelectedVertices; }
            set
            {
                Controller.DrawSelectedVertices = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Inverts the selection of vertices <see cref="GridEditorState.SelectedVertices"/>
        /// </summary>
        public bool InvertVerticesSelection
        {
            get { return Controller.InvertVerticesSelection; }
            set
            {
                Controller.InvertVerticesSelection = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Size of the grid vertices (when in edit mode <see cref="IsEditing"/>)
        /// </summary>
        public float VertexSize
        {
            get { return Controller.VertexSize; }
            set
            {
                Controller.VertexSize = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Command for executing the <see cref="SelectedGenerateGridWorkFlowType"/>
        /// </summary>
        public ICommand GenerateUsingCurrentWorkFlowCommand { get; private set; }

        /// <summary>
        /// Command for starting a map tool (using <see cref="MapToolType"/> command parameter
        /// </summary>
        public ICommand SelectMapToolCommand { get; private set; }

        /// <summary>
        /// Command for starting polygon change tool
        /// </summary>
        public ICommand ChangePolygonCommand { get; private set; }

        /// <summary>
        /// Command for refining a polygon
        /// </summary>
        public ICommand RefinePolygonCommand { get; private set; }

        /// <summary>
        /// Command for offsetting a polygon
        /// </summary>
        public ICommand OffsetPolygonCommand { get; private set; }

        /// <summary>
        /// Deletes the selected vertices
        /// </summary>
        public ICommand DeleteVerticesCommand { get; private set; }

        /// <summary>
        /// Deletes mesh with options
        /// </summary>
        public ICommand DeleteMeshWithOptionsCommand { get; private set; }
        

        public ICommand FlipEdgesCommand { get; private set; }

        public ICommand MergeTwoVerticesCommand { get; private set; }

        public ICommand MergeVerticesCommand { get; private set; }

        public ICommand OrthogonalizeCommand { get; private set; }

        public ICommand MakeGridCommand { get; private set; }

        public ICommand MakeGridFromSplinesCommand { get; private set; }

        public ICommand MakeGridFromSplinesOrthogonalCommand { get; private set; }

        public ICommand MakeTriangularGridInPolygon { get; private set; }

        public ICommand MakeTriangularGridFromSamples { get; private set; }

        public ICommand RefineGridBasedOnSampleCommand { get; private set; }

        public ICommand RefineGridBasedOnPolygonCommand { get; private set; }

        public ICommand GetMeshBoundaryPolygon { get; private set; }

        public ICommand ShowOrthogonality { get; private set; }

        public ICommand ShowSmoothness { get; private set; }

        /// <summary>
        /// Starts the import of the provided <see cref="ImportExportType"/>
        /// </summary>
        public ICommand ImportCommand { get; private set; }

        /// <summary>
        /// Starts the import of the provided <see cref="ImportExportType"/>
        /// </summary>
        public ICommand ExportCommand { get; private set; }

        /// <summary>
        /// Switches the current workflow
        /// </summary>
        public ICommand SwitchWorkFlowCommand { get; private set; }

        /// <summary>
        /// Is the mask tool currently active
        /// </summary>
        public bool PolygonToolEnabled
        {
            get { return Controller?.SelectedMapToolType == MapToolType.Polygon; }
        }

        public bool ChangePolygonToolEnabled
        {
            get { return Controller?.SelectedMapToolType == MapToolType.ChangePolygon; }
        }

        /// <summary>
        /// Is the spline tool currently active
        /// </summary>
        public bool SplineToolEnabled
        {
            get { return Controller?.SelectedMapToolType == MapToolType.Spline; }
        }

        /// <summary>
        /// Is the spline tool currently active
        /// </summary>
        public bool LandBoundariesToolEnabled
        {
            get { return Controller?.SelectedMapToolType == MapToolType.LandBoundaries; }
        }

        /// <summary>
        /// Is the insert edges tools enabled
        /// </summary>
        public bool InsertEdgesToolEnabled
        {
            get { return Controller?.SelectedMapToolType == MapToolType.InsertEdges; }
        }

        public OrthogonalizationParameters OrthogonalizationParameters
        {
            get { return Controller.OrthogonalizationParameters; }
        }

        public MakeGridParameters MakeGridParameters
        {
            get { return Controller.MakeGridParameters; }
        }

        public CurvilinearParameters CurvilinearParameters
        {
            get { return Controller.CurvilinearParameters; }
        }

        public SplinesToCurvilinearParameters SplinesToCurvilinearParameters
        {
            get { return Controller.SplinesToCurvilinearParameters; }
        }

        public InterpolationParameters InterpolationParameters
        {
            get { return Controller.InterpolationParameters; }
        }

        public SamplesRefineParameters SamplesRefineParameters
        {
            get { return Controller.SamplesRefineParameters; }
        }

        public double PolygonRefinementDistance
        {
            get => Controller.PolygonRefinementDistance;
            set
            {
                Controller.PolygonRefinementDistance = value;
                OnPropertyChanged();
            }
        }

        public double PolygonOffsetDistance
        {
            get => Controller.PolygonOffsetDistance;
            set
            {
                Controller.PolygonOffsetDistance = value;
                OnPropertyChanged();
            }
        }

        public bool IsTriangulationRequired
        {
            get { return Controller.IsTriangulationRequired; }
            set
            {
                Controller.IsTriangulationRequired = value;
                OnPropertyChanged();
            }
        }

        public bool IsAccountingForLandBoundariesRequired
        {
            get { return Controller.IsAccountingForLandBoundariesRequired; }
            set
            {
                Controller.IsAccountingForLandBoundariesRequired = value;
                OnPropertyChanged();
            }
        }

        public string CurrentProgressText
        {
            get { return currentProgressText; }
            set
            {
                currentProgressText = value;
                OnPropertyChanged();
            }
        }

        public float CurrentProgress
        {
            get { return currentProgress; }
            set
            {
                currentProgress = value;
                OnPropertyChanged();
            }
        }

        [TypeConverter(typeof(EnumDescriptionAttributeTypeConverter))]
        public ProjectToLandBoundaryOptions ProjectToLandBoundaryOption
        {
            get { return Controller.ProjectToLandBoundaryOption; }
            set
            {
                Controller.ProjectToLandBoundaryOption = value;
                OnPropertyChanged();
            }
        }

        [TypeConverter(typeof(EnumDescriptionAttributeTypeConverter))]
        public DeleteMeshOptions DeleteMeshOption
        {
            get { return Controller.DeleteMeshOption; }
            set
            {
                Controller.DeleteMeshOption = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Controller containing grid editing logic 
        /// </summary>
        internal IGridEditorController Controller { get; }

        /// <summary>
        /// Gets the file path for the selected <see cref="ImportExportAction"/>
        /// </summary>
        internal Func<string, ImportExportAction, string> GetFilePath { private get; set; }

        /// <inheritdoc/>
        public void RefreshState()
        {
            OnPropertyChanged(nameof(PolygonToolEnabled));
            OnPropertyChanged(nameof(ChangePolygonToolEnabled));
            OnPropertyChanged(nameof(SplineToolEnabled));
            OnPropertyChanged(nameof(LandBoundariesToolEnabled));
            OnPropertyChanged(nameof(InsertEdgesToolEnabled));

            OnPropertyChanged(nameof(ProjectToLandBoundaryOption));
            OnPropertyChanged(nameof(DeleteMeshOption));
            OnPropertyChanged(nameof(IsAccountingForLandBoundariesRequired));
            OnPropertyChanged(nameof(IsTriangulationRequired));

            OnPropertyChanged(nameof(SplinesToCurvilinearParameters));
            OnPropertyChanged(nameof(CurvilinearParameters));
            OnPropertyChanged(nameof(PolygonRefinementDistance));
            OnPropertyChanged(nameof(MakeGridParameters));
            OnPropertyChanged(nameof(OrthogonalizationParameters));

            OnPropertyChanged(nameof(VertexSize));
            OnPropertyChanged(nameof(DrawSelectedVertices));
            OnPropertyChanged(nameof(DrawFilledInPolygons));
            OnPropertyChanged(nameof(ShowVertices));

            OnPropertyChanged(nameof(IsEditing));
        }

        public bool IsTaskRunning
        {
            get { return isTaskRunning; }
            set
            {
                isTaskRunning = value;
                OnPropertyChanged();
                OnPropertyChanged(nameof(CanDoActions));
            }
        }

        public void UpdateProgress(float percentage, string progressText)
        {
            CurrentProgress = percentage;
            CurrentProgressText = progressText;
        }

        /// <inheritdoc/>
        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        private void SetCommands()
        {
            GenerateUsingCurrentWorkFlowCommand = new RelayCommand(o => Controller.GenerateGridForSelectedWorkflow(), o => IsEditing);

            SelectMapToolCommand = new RelayCommand(o => { Controller.EnableMapToolByType((MapToolType) o); });
            ChangePolygonCommand = new RelayCommand(o =>
            {
                Controller.EnableMapToolByType(MapToolType.ChangePolygon);
            });
            RefinePolygonCommand = new RelayCommand(o => { Controller.RefineSelectedPolygons(); });

            OffsetPolygonCommand = new RelayCommand(o => { Controller.AddOffsettedPolygon(); });

            FlipEdgesCommand = new RelayCommand(o => { Controller.FlipEdges(); }, o => IsEditing);

            DeleteVerticesCommand = new RelayCommand(o => { Controller.DeleteSelectedVertices(); }, o => IsEditing);

            DeleteMeshWithOptionsCommand = new RelayCommand(o => { Controller.DeleteMeshWithOptions(); }, o => IsEditing);
            MergeTwoVerticesCommand = new RelayCommand(o => { Controller.MergeTwoVertices(); },
                o => IsEditing && Controller.HasSelectedVertices);
            MergeVerticesCommand = new RelayCommand(o => { Controller.MergeVertices(); },
                o => IsEditing && Controller.HasSelectedVertices);

            OrthogonalizeCommand = new RelayCommand(o => { Controller.Orthogonalize(); }, o => IsEditing);

            MakeGridCommand = new RelayCommand(o => { Controller.MakeGrid(); }, o => IsEditing);

            MakeGridFromSplinesCommand = new RelayCommand(o => { Controller.MakeGridFromSplines(); }, o => IsEditing);
            MakeGridFromSplinesOrthogonalCommand =
                new RelayCommand(o => { Controller.MakeGridFromSplinesOrthogonal(); }, o => IsEditing);

            MakeTriangularGridInPolygon =
                new RelayCommand(o => { Controller.MakeTriangularGridInPolygon(); }, o => IsEditing);


            MakeTriangularGridFromSamples=
                new RelayCommand(o => { Controller.MakeTriangularGridFromSamples(); }, o => IsEditing);

            RefineGridBasedOnSampleCommand = new RelayCommand(o => { Controller.RefineGridBasedOnSamples(); },
                o => IsEditing && Controller.HasSamples);

            RefineGridBasedOnPolygonCommand =
                new RelayCommand(o => { Controller.RefineGridBasedOnPolygon(); }, o => IsEditing);

            GetMeshBoundaryPolygon = new RelayCommand(o => { Controller.GetMeshBoundaryPolygon(); }, o => IsEditing);


            ShowOrthogonality = new RelayCommand(o => { Controller.ShowOrthogonality(); }, o => IsEditing);

            ShowSmoothness = new RelayCommand(o => { Controller.ShowSmoothness(); }, o => IsEditing);

            ImportCommand = new RelayCommand(o =>
            {
                SelectedImportType = (ImportExportType)o;
                ImportExportGridStateData(ImportExportAction.Import, SelectedImportType);
            });

            ExportCommand = new RelayCommand(o =>
            {
                SelectedExportType = (ImportExportType)o;
                ImportExportGridStateData(ImportExportAction.Export, SelectedExportType);
            });

            SwitchWorkFlowCommand = new RelayCommand(o => { SelectedGenerateGridWorkFlowType = (GenerateGridWorkFlowType) o; });

            CancelTaskCommand = new RelayCommand(o => Controller.RequestCancelation(), o => IsTaskRunning);
        }

        public ICommand CancelTaskCommand { get; private set; }

        private string GetFileFilter(ImportExportType importExportType)
        {
            switch (importExportType)
            {
                case ImportExportType.Grid:
                    return "Net file|*.nc";
                case ImportExportType.Polygons:
                    return "Shape file|*.shp|pol file|*.pol";
                case ImportExportType.Spline:
                    return "Shape file|*.shp|spl file|*.spl";
                case ImportExportType.LandBoundary:
                    return "Shape file|*.shp";
                case ImportExportType.Samples:
                    return "xyz file|*.xyz|ascii file|*.asc";
                default:
                    throw new ArgumentOutOfRangeException(nameof(importExportType), importExportType, null);
            }
        }

        private void ImportExportGridStateData(object o, ImportExportType importExportType)
        {
            if (!(o is ImportExportAction action)) return;

            var filePath = GetFilePath?.Invoke(GetFileFilter(importExportType), action);
            if (filePath == null) return;

            Controller.ImportExport(action, importExportType, filePath);
        }

        public bool IsPolygonToolVisible
        {
            get { return SelectedGenerateGridWorkFlowType == GenerateGridWorkFlowType.RegularWithPolygon ||
                         SelectedGenerateGridWorkFlowType == GenerateGridWorkFlowType.TriangularWithPolygon; }
        }

        public bool IsSplineToolVisible
        {
            get { return SelectedGenerateGridWorkFlowType == GenerateGridWorkFlowType.CurveLinearOrthogonal ||
                         SelectedGenerateGridWorkFlowType == GenerateGridWorkFlowType.CurveLinearTransfinite;
            }
        }
    }
}