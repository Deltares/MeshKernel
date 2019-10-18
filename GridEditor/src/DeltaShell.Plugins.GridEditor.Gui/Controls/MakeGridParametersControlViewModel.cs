using System;
using System.ComponentModel;
using System.IO;
using System.Runtime.CompilerServices;
using System.Text;
using System.Windows.Input;
using DelftTools.Controls.Wpf.Commands;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.SharpMapGis;
using log4net;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    public class MakeGridParametersControlViewModel : INotifyPropertyChanged
    {
        private static readonly ILog Log = LogManager.GetLogger(typeof(MakeGridParametersControlViewModel));

        private MakeGridParameters gridParameters;
        private string asciiFilePath;
        private bool forcePointsOnCellCenters;
        
        private double gridBlockSize;
        private double xGridBlockSize;
        private double yGridBlockSize;
        private double originXCoordinate;
        private double originYCoordinate;
        private int numberOfRows;
        private int numberOfColumns;

        public MakeGridParametersControlViewModel()
        {
            ReadAsciiFileDimensionsCommand = new RelayCommand(o =>
            {
                AsciiFilePath = GetFilePath?.Invoke();
                ReadParameters();
            });

            SetDimensionsFromAsciiCommand = new RelayCommand(o => CopyStateToParameters(), o => AscFileExists);
        }

        /// <summary>
        /// <see cref="GridParameters"/> to edit
        /// </summary>
        public MakeGridParameters GridParameters
        {
            get { return gridParameters; }
            set
            {
                gridParameters = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Path ot the file used for making a grid
        /// </summary>
        public string AsciiFilePath
        {
            get { return asciiFilePath; }
            set
            {
                asciiFilePath = value;
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// Force the points of the ascii file on the grid centers
        /// </summary>
        public bool ForcePointsOnCellCenters
        {
            get { return forcePointsOnCellCenters; }
            set
            {
                forcePointsOnCellCenters = value;
                FileParametersString = GetParametersString();
                OnPropertyChanged();
                OnPropertyChanged(nameof(FileParametersString));
            }
        }

        internal Func<string> GetFilePath { private get; set; }

        /// <summary>
        /// Command to set the specified <see cref="AsciiFilePath"/> and read the dimensions
        /// </summary>
        public ICommand ReadAsciiFileDimensionsCommand { get; }

        /// <summary>
        /// Sets the read dimensions to the <see cref="GridParameters"/>
        /// </summary>
        public ICommand SetDimensionsFromAsciiCommand { get; }

        /// <summary>
        /// String showing the currently read dimensions
        /// </summary>
        public string FileParametersString { get; private set; } = "";

        private bool AscFileExists
        {
            get { return !string.IsNullOrEmpty(AsciiFilePath) && File.Exists(AsciiFilePath); }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        private void OnPropertyChanged([CallerMemberName] string propertyName = null)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(propertyName));
        }

        private void ReadParameters()
        {
            if (string.IsNullOrEmpty(AsciiFilePath) || !File.Exists(AsciiFilePath))
            {
                Log.Error($"Can not find file {AsciiFilePath}");
                return;
            }

            ReadDimensionsFromAscii();
            FileParametersString = GetParametersString();
            OnPropertyChanged(nameof(FileParametersString));
        }

        private void ReadDimensionsFromAscii()
        {
            var grid = RegularGridCoverageLoader.LoadFromFile(AsciiFilePath);

            gridBlockSize = grid.DeltaX;
            xGridBlockSize = grid.DeltaX;
            yGridBlockSize = grid.DeltaX;

            originXCoordinate = grid.Origin.X;
            originYCoordinate = grid.Origin.Y;

            numberOfRows = grid.SizeY;
            numberOfColumns = grid.SizeX;
        }

        private string GetParametersString()
        {
            var stringBuilder = new StringBuilder();

            stringBuilder.AppendLine($"Type = {GridTypeOptions.Square}");
            stringBuilder.AppendLine($"BlockSize = {gridBlockSize}");
            stringBuilder.AppendLine($"X BlockSize = {xGridBlockSize}");
            stringBuilder.AppendLine($"Y BlockSize = {yGridBlockSize}");
            stringBuilder.AppendLine($"Origin X Coordinate = {GetOriginXCoordinate()}");
            stringBuilder.AppendLine($"Origin Y Coordinate = {GetOriginYCoordinate()}");
            stringBuilder.AppendLine($"Number Of Rows  = {GetNumberOfRows()}");
            stringBuilder.AppendLine($"Number Of Columns = {GetNumberOfColumns()}");

            return stringBuilder.ToString();
        }

        private void CopyStateToParameters()
        {
            GridParameters.GridBlockSize = gridBlockSize;
            GridParameters.XGridBlockSize = xGridBlockSize;
            GridParameters.YGridBlockSize = yGridBlockSize;

            GridParameters.OriginXCoordinate = GetOriginXCoordinate();
            GridParameters.OriginYCoordinate = GetOriginYCoordinate();
            GridParameters.NumberOfRows = GetNumberOfRows();
            GridParameters.NumberOfColumns = GetNumberOfColumns();
        }

        private int GetNumberOfColumns()
        {
            return ForcePointsOnCellCenters
                ? numberOfColumns
                : numberOfColumns - 1;
        }

        private double GetOriginXCoordinate()
        {
            return ForcePointsOnCellCenters
                ?  originXCoordinate - xGridBlockSize / 2.0
                : originXCoordinate;
        }

        private double GetOriginYCoordinate()
        {
            return ForcePointsOnCellCenters
                ? originYCoordinate - yGridBlockSize/ 2.0
                : originYCoordinate;
        }

        private int GetNumberOfRows()
        {
            return ForcePointsOnCellCenters
                ? numberOfRows
                : numberOfRows - 1;
        }
    }
}