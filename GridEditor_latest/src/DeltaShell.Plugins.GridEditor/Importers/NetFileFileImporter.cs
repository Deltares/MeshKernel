using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Drawing;
using DelftTools.Shell.Core;
using DelftTools.Utils.Collections.Extensions;
using DeltaShell.Plugins.SharpMapGis.ImportExport;
using NetTopologySuite.Extensions.Grids;

namespace DeltaShell.Plugins.GridEditor.Importers
{
    public class NetFileFileImporter : IFileImporter
    {
        public string Name
        {
            get { return "Net file"; }
        }

        public string Category
        {
            get { return "Grid importers"; }
        }

        public string Description
        {
            get { return "Net file"; }
        }

        public Bitmap Image { get; }

        public IEnumerable<Type> SupportedItemTypes
        {
            get { yield return typeof(UnstructuredGrid); }
        }

        public bool CanImportOnRootLevel
        {
            get { return true; }
        }

        public string FileFilter
        {
            get { return "Net file|*.nc"; }
        }

        [ExcludeFromCodeCoverage]
        public string TargetDataDirectory { get; set; }

        [ExcludeFromCodeCoverage]
        public bool ShouldCancel { get; set; }

        [ExcludeFromCodeCoverage]
        public ImportProgressChangedDelegate ProgressChanged { get; set; }

        public bool OpenViewAfterImport
        {
            get { return true; }
        }

        public bool CanImportOn(object targetObject)
        {
            return targetObject is Project ||
                   targetObject is UnstructuredGrid;
        }

        public object ImportItem(string path, object target = null)
        {
            if (target is Project)
            {
                return ImportGridFromFile(path);
            }

            if (target is UnstructuredGrid unstructuredGrid)
            {
                unstructuredGrid.Clear();

                var importGrid = ImportGridFromFile(path);
                unstructuredGrid.Vertices.AddRange(importGrid.Vertices);
                unstructuredGrid.Edges.AddRange(importGrid.Edges);

                return importGrid;
            }

            return null;
        }

        internal Func<string, UnstructuredGrid> ImportGridFromFile { get; set; } = path =>
        {
            return NetFileImporter.ImportGrid(path);
        };
    }
}