using System.ComponentModel;
using System.IO;
using System.Runtime.CompilerServices;
using System.Text;
using DeltaShell.Plugins.GridEditor.GridGeomStateful.Api;
using DeltaShell.Plugins.SharpMapGis;
using GeoAPI.Extensions.Coverages;
using log4net;

namespace DeltaShell.Plugins.GridEditor.Gui.Controls
{
    public class MakeGridParametersControlViewModel : INotifyPropertyChanged
    {
        private static readonly ILog Log = LogManager.GetLogger(typeof(MakeGridParametersControlViewModel));

        private MakeGridParameters gridParameters;
        private string asciiFilePath;

        private double gridBlockSize;
        private double xGridBlockSize;
        private double yGridBlockSize;
        private double originXCoordinate;
        private double originYCoordinate;
        private int numberOfRows;
        private int numberOfColumns;
        private bool forcePointsOnCellCenters = true;
        private bool isBasedOnPolygon;
        
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
                asciiFilePath = value?.ToLower();
                ReadParameters();
                OnPropertyChanged();
            }
        }

        /// <summary>
        /// String showing the currently read dimensions
        /// </summary>
        public string FileParametersString { get; private set; } = "";

        public bool IsBasedOnPolygon
        {
            get { return isBasedOnPolygon; }
            set
            {
                isBasedOnPolygon = value;
                OnPropertyChanged();
            }
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
            forcePointsOnCellCenters = grid.PointPosition == RegularGridPointPosition.OnCellCentre;

            originXCoordinate = grid.Origin.X;
            originYCoordinate = grid.Origin.Y;

            numberOfRows = grid.SizeY;
            numberOfColumns = grid.SizeX;

            CopyStateToParameters();
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
            return numberOfColumns;
        }

        private double GetOriginXCoordinate()
        {
            return forcePointsOnCellCenters
                ?  originXCoordinate - xGridBlockSize / 2.0
                : originXCoordinate;
        }

        private double GetOriginYCoordinate()
        {
            return forcePointsOnCellCenters
                ? originYCoordinate - yGridBlockSize/ 2.0
                : originYCoordinate;
        }

        private int GetNumberOfRows()
        {
            return numberOfRows;

        }
    }
}